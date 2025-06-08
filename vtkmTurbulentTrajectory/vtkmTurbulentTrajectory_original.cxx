#include <cstdlib>
#include <fstream>
#include <highfive/highfive.hpp>
#include <iostream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <vtkm/CellShape.h> 
#include <vtkm/Math.h>
#include <vtkm/Matrix.h>
#include <vtkm/Particle.h>
#include <vtkm/Types.h>
#include <vtkm/VecFromPortal.h>
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleRandomUniformReal.h>
#include <vtkm/cont/CellLocatorGeneral.h>
#include <vtkm/cont/CellSetExplicit.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/filter/resampling/Probe.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/random/Philox.h>
#include <vtkm/worklet/WorkletMapField.h>

class LangevinWorklet : public vtkm::worklet::WorkletMapField {
public:
  LangevinWorklet(const vtkm::Float32 &stepSize, const vtkm::Vec3f &minCoords,
                  const vtkm::Vec3f &maxCoords)
      : stepSize_(stepSize), minCoords_(minCoords), maxCoords_(maxCoords) {}

  using ControlSignature = void(FieldIn, FieldIn, FieldIn, FieldIn, FieldInOut);

  VTKM_EXEC void operator()(const vtkm::Vec3f &U,
                            const vtkm::Vec<vtkm::FloatDefault, 6> &RST,
                            const vtkm::Vec3f &Normal3,
                            const vtkm::Vec3f &randomPositionAtEntrance,
                            vtkm::Vec3f &p0) const {

    vtkm::Vec3f normalizedNormal3 = vtkm::Normal(Normal3);

    vtkm::Matrix<vtkm::FloatDefault, 3, 3> RST_;

    RST_(0, 0) = RST[0];
    RST_(0, 1) = RST[3];
    RST_(0, 2) = RST[5];
    RST_(1, 0) = RST[3];
    RST_(1, 1) = RST[1];
    RST_(1, 2) = RST[4];
    RST_(2, 0) = RST[5];
    RST_(2, 1) = RST[3];
    RST_(2, 2) = RST[2];

    vtkm::Vec3f vt = vtkm::MatrixMultiply(RST_, normalizedNormal3);

    p0 = p0 + (U + vt * vtkm::Sqrt(1 - vtkm::Exp(-2 * stepSize_))) * stepSize_;

    if (!(p0[0] >= minCoords_[0] && p0[0] <= maxCoords_[0] &&
          p0[1] >= minCoords_[1] && p0[1] <= maxCoords_[1] &&
          p0[2] >= minCoords_[2] && p0[2] <= maxCoords_[2])) {
      p0 = randomPositionAtEntrance;
    }
  }

private:
  const vtkm::Float32 stepSize_;
  const vtkm::Vec3f minCoords_;
  const vtkm::Vec3f maxCoords_;
};

int main(int argc, char **argv) {

  // VTKm initialization (any device)
  auto opts = vtkm::cont::InitializeOptions::DefaultAnyDevice;
  auto config = vtkm::cont::Initialize(argc, argv, opts);

  std::string settings_file = std::string(argv[1]);

  std::ifstream settingsFile(settings_file);
  if (!settingsFile.is_open()) {
    std::cerr << "Failed to open settings file." << std::endl;
    return 1;
  }

  nlohmann::json settings;
  settingsFile >> settings;

  std::string dataFile = settings["dataFile"];
  std::string fieldName1 = settings["fieldNameU"];
  std::string fieldName2 = settings["fieldNameRST"];
  int numStepsToWrite = settings["numStepsToWrite"];
  std::string outputFolder = settings["outputFolder"];
  std::string seedName = settings["seedName"];
  std::string entranceSeedName = settings["entranceName"];
  std::string h5filename = settings["h5filename"];
  vtkm::Float32 stepSize = settings["stepSize"];
  vtkm::Float32 delta_t = settings["delta_t"];

  HighFive::File file(outputFolder + "/" + h5filename,
                      HighFive::File::ReadWrite | HighFive::File::Create |
                          HighFive::File::Truncate);

  if (argc < 1) {
    std::cerr << "Usage:\n"
              << "vtkmTurbulentTrajectory settings_file.json\n"
              << std::endl;
    std::cerr << "Example:\n";
    std::cerr << "./vtkmTurbulentTrajectory advection_settings.json\n";
    std::cerr << "where vtk-m options are: " << std::endl
              << config.Usage << std::endl;
    exit(EXIT_FAILURE);
  }

  // Read seeds from a file and store them in a vtkm::cont::ArrayHandle
  std::ifstream inputFile(seedName);
  if (!inputFile.is_open()) {
    std::cerr << "Failed to open coordinates file." << std::endl;
    return 1;
  }

  std::vector<vtkm::Vec3f> coordinates;
  std::string line;
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    vtkm::Vec3f coordinate;
    if (!(iss >> coordinate[0] >> coordinate[1] >> coordinate[2])) {
      std::cerr << "Failed to read coordinate from file." << std::endl;
      return 1;
    }
    coordinates.push_back(coordinate);
  }

  // Print the max and min values of the domain of the seeds
  vtkm::Vec3f minCoords =
      vtkm::Vec3f(std::numeric_limits<vtkm::Float32>::max());
  vtkm::Vec3f maxCoords =
      vtkm::Vec3f(std::numeric_limits<vtkm::Float32>::lowest());
  for (const auto &coordinate : coordinates) {
    minCoords[0] = std::min(minCoords[0], coordinate[0]);
    minCoords[1] = std::min(minCoords[1], coordinate[1]);
    minCoords[2] = std::min(minCoords[2], coordinate[2]);
    maxCoords[0] = std::max(maxCoords[0], coordinate[0]);
    maxCoords[1] = std::max(maxCoords[1], coordinate[1]);
    maxCoords[2] = std::max(maxCoords[2], coordinate[2]);
  }
  std::cout << "Minimum coordinates (model)    : " << minCoords << std::endl;
  std::cout << "Maximum coordinates (model)    : " << maxCoords << std::endl;

  // Read reseed location from a file and store them in a
  // vtkm::cont::ArrayHandle
  std::ifstream entranceSeedsFile(entranceSeedName);
  if (!entranceSeedsFile.is_open()) {
    std::cerr << "Failed to open coordinates file." << std::endl;
    return 1;
  }

  std::vector<vtkm::Vec3f> coordinatesEntrance;
  std::string lineEntrance;
  while (std::getline(entranceSeedsFile, lineEntrance)) {
    std::istringstream iss(lineEntrance);
    vtkm::Vec3f coordinate;
    if (!(iss >> coordinate[0] >> coordinate[1] >> coordinate[2])) {
      std::cerr << "Failed to read coordinate from entrance file." << std::endl;
      return 1;
    }
    coordinatesEntrance.push_back(coordinate);
  }

  vtkm::cont::ArrayHandle<vtkm::Vec3f> coordinatesEntranceArray;
  coordinatesEntranceArray.Allocate(
      static_cast<vtkm::Id>(coordinatesEntrance.size()));
  auto portal = coordinatesEntranceArray.WritePortal();
  for (vtkm::Id i = 0; i < coordinatesEntranceArray.GetNumberOfValues(); ++i) {
    portal.Set(i, coordinatesEntrance[i]);
  }

  vtkm::cont::Timer Timer;
  Timer.Start();
  std::cout << "Reading data file... " << std::endl;
  vtkm::io::VTKDataSetReader reader(dataFile);
  vtkm::cont::DataSet ds0 = reader.ReadDataSet();

  Timer.Stop();
  std::cout << "Data file read time: " << Timer.GetElapsedTime() << " seconds"
            << std::endl;
  Timer.Reset();

  std::vector<vtkm::cont::ArrayHandle<vtkm::Vec3f>> updatedPositionsVector;

  int numSteps = static_cast<int>(std::abs(delta_t / stepSize));

  // Build the probe positions for the first timestep and store them in a file
  vtkm::cont::DataSetBuilderExplicitIterative dataSetBuilder0;
  for (std::vector<vtkm::Vec3f>::size_type i = 0; i < coordinates.size(); i++) {
    dataSetBuilder0.AddPoint(coordinates[i]);
    dataSetBuilder0.AddCell(vtkm::CELL_SHAPE_VERTEX);
    dataSetBuilder0.AddCellPoint(i);
  }
  vtkm::cont::DataSet initialPositions = dataSetBuilder0.Create();

  std::cout << "Initial positions created." << std::endl;

  // Build the probeFilter
  vtkm::filter::resampling::Probe probeFilter;
  probeFilter.SetFieldsToPass({fieldName1, fieldName2});

  // Loop over timesteps (some particles may leave the domain and need to be
  // reseeded)
  vtkm::cont::DataSet outputDataset = initialPositions;

  std::cout << "Starting loop over timesteps..." << std::endl;

  Timer.Start();
  for (int i = 0; i <= numSteps; i++) {
    vtkm::cont::ArrayHandle<vtkm::Vec3f> positions;
    vtkm::cont::ArrayCopy(outputDataset.GetCoordinateSystem().GetData(),
                          positions);

    // Create vector of positions to sample Velocity field and RST field
    vtkm::Id numPts = positions.GetNumberOfValues();


    // Probe velocity field at the current position
    probeFilter.SetGeometry(outputDataset);
    probeFilter.SetInvalidValue(0.0f);
    vtkm::cont::DataSet probedData = probeFilter.Execute(ds0);

    // std::cout << "Probed data size: " << probedData.GetNumberOfCells()
    //           << std::endl;

    // Get the field data from probedData at every position
    // std::cout << "Getting velocity field data... " << std::endl;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> velocityField;
    probedData.GetField("U").GetData().AsArrayHandle(velocityField);

    // Random normal vector using ArrayHandleRandomUniformReal
    vtkm::Id NumPoints = probedData.GetNumberOfCells();

    auto Normal3 = vtkm::cont::make_ArrayHandleCompositeVector(
        vtkm::cont::ArrayHandleRandomUniformReal<vtkm::FloatDefault>(NumPoints),
        vtkm::cont::ArrayHandleRandomUniformReal<vtkm::FloatDefault>(NumPoints),
        vtkm::cont::ArrayHandleRandomUniformReal<vtkm::FloatDefault>(
            NumPoints));

    // Random position at entrance
    std::vector<vtkm::Vec3f> randomPositionsAtEntrance_;
    for (vtkm::Id k = 0; k < NumPoints; k++) {
      randomPositionsAtEntrance_.push_back(
          coordinatesEntrance[rand() % coordinatesEntrance.size()]);
    }
    vtkm::cont::ArrayHandle<vtkm::Vec3f> randomPositionsAtEntrance =
        vtkm::cont::make_ArrayHandle(randomPositionsAtEntrance_, vtkm::CopyFlag::On);

    // std::cout << "Getting RST field data... " << std::endl;
    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 6>> RST;
    probedData.GetField("RST").GetData().AsArrayHandle(RST);


    vtkm::cont::Invoker invoke;
    auto LangevinStep = LangevinWorklet(stepSize, minCoords, maxCoords);
    invoke(LangevinStep, velocityField, RST, Normal3, randomPositionsAtEntrance, positions);

    vtkm::cont::DataSetBuilderExplicitIterative dataSetBuilder;
    auto updatedPositionPortal = positions.ReadPortal();
    std::vector<std::vector<double>> outputVector;

    for (vtkm::Id k = 0; k < numPts; k++) {
      vtkm::Vec3f updatedPosition = updatedPositionPortal.Get(k);
      dataSetBuilder.AddPoint(updatedPosition);
      dataSetBuilder.AddCell(vtkm::CELL_SHAPE_VERTEX);
      dataSetBuilder.AddCellPoint(k);
      outputVector.push_back(
          {updatedPosition[0], updatedPosition[1], updatedPosition[2]});
    }
    outputDataset = dataSetBuilder.Create();

    // Code to execute if numSteps is a multiple of numStepsToWrite
    if (i % numStepsToWrite == 0 && i != 0) {
      std::cout << "Writing output for iteration " << i << std::endl;
      std::cout << "Number of particles: " << numPts << std::endl;
      file.createDataSet("Iteration_" + std::to_string(i), outputVector);
      std::string outputTime = outputFolder + "/output_" + std::to_string(i);
      vtkm::io::VTKDataSetWriter writer(outputTime + ".vtk");
      writer.WriteDataSet(outputDataset);
    }
    outputVector.clear();
  }
  Timer.Stop();
  std::cout << "Total simulation time: " << Timer.GetElapsedTime() << " seconds"
            << std::endl;

  return 0;
}
