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

  using ControlSignature = void(FieldIn, FieldIn, FieldIn, FieldIn, FieldInOut, FieldInOut);

  VTKM_EXEC void operator()(const vtkm::Vec3f &U,
                            const vtkm::Vec<vtkm::FloatDefault, 6> &RST,
                            const vtkm::Vec3f &Normal3,
                            const vtkm::Vec3f &randomPositionAtEntrance,
                            vtkm::Vec3f &p0,
                            vtkm::UInt8 &activeFlag) const {

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

    p0 = p0 + (U + vt * vtkm::Sqrt(1 - vtkm::Exp(-2 * stepSize_))) * stepSize_; // Turbulent
    // p0 = p0 + U * stepSize_; // Non-turbulent

    if (!(p0[0] >= minCoords_[0] && p0[0] <= maxCoords_[0] &&
          p0[1] >= minCoords_[1] && p0[1] <= maxCoords_[1] &&
          p0[2] >= minCoords_[2] && p0[2] <= maxCoords_[2])) {
      p0 = randomPositionAtEntrance;
      activeFlag = 1; // Mark as reseeded
    }
  }

private:
  const vtkm::Float32 stepSize_;
  const vtkm::Vec3f minCoords_;
  const vtkm::Vec3f maxCoords_;
};

// NEW ----------->
float average(std::vector<float> const& v){
    if(v.empty()){
        return 0;
    }
    auto const count = static_cast<float>(v.size());
    return std::reduce(v.begin(), v.end()) / count;
}
void subtractMean(std::vector<float>& v) {
    float mean = average(v);
    for (auto& element : v) {
        element -= mean;
    }
}
// <---------- NEW 

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
  bool turbulent = settings["turbulent"];

  // NEW ----------->
  std::string phantomFileName = settings["phantomFileName"];
  float T1 = settings["T1"];
  float T2 = settings["T2"];
  float PD = settings["PD"];
  // <---------- NEW 

  // HighFive::File file(outputFolder + "/" + h5filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

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

  // NEW ----------->     
  int Nspins = coordinates.size();
  std::cout << "Number of spins: " << Nspins << std::endl;

  std::vector<float> xVector;
  std::vector<float> yVector;
  std::vector<float> zVector;
  std::vector<float> t1Vector;
  std::vector<float> t2Vector;
  std::vector<float> pdVector;
  // dx, dy, dz and spin_reset matrices with Nspins filas:
  std::vector<std::vector<float>> dx;
  std::vector<std::vector<float>> dy;
  std::vector<std::vector<float>> dz;
  std::vector<std::vector<int>> spin_reset;

  std::vector<float> t      = {0.0, delta_t};
  std::vector<float> t_unit = {0.0, 1.0};
  // <---------- NEW 

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

    // NEW ----------->    
    xVector.push_back(coordinates[i][0]);
    yVector.push_back(coordinates[i][1]);
    zVector.push_back(coordinates[i][2]);
    t1Vector.push_back(T1);
    t2Vector.push_back(T2);
    pdVector.push_back(PD);
    // <---------- NEW 
  }
  vtkm::cont::DataSet initialPositions = dataSetBuilder0.Create();

  std::cout << "Initial positions created." << std::endl;

  // Build the probeFilter
  vtkm::filter::resampling::Probe probeFilter;
  probeFilter.SetFieldsToPass({fieldName1, fieldName2});

  // NEW ----------->
  std::vector<float> rowX(Nspins, 0.0f);
  std::vector<float> rowY(Nspins, 0.0f);
  std::vector<float> rowZ(Nspins, 0.0f);
  std::vector<int> rowReset(Nspins, 0);
  dx.push_back(rowX);
  dy.push_back(rowY);
  dz.push_back(rowZ);
  spin_reset.push_back(rowReset);

  vtkm::cont::ArrayHandle<vtkm::UInt8> flags; // flags for every particle to store the previous "inside" boolean. Initialize to an array of zeros
  flags.Allocate(Nspins); flags.Fill(0);      // Initialize all flags to 0
  // <---------- NEW 

  // Loop over timesteps (some particles may leave the domain and need to be reseeded)
  vtkm::cont::DataSet outputDataset = initialPositions;

  std::cout << "Starting loop over timesteps..." << std::endl;

  if (turbulent) {
      std::cout << "Turbulent... " << std::endl;
  } else {
      std::cout << "Non-Turbulent... " << std::endl;
  }

  Timer.Start();
  for (int i = 0; i <= numSteps; i++) {
    vtkm::cont::ArrayHandle<vtkm::Vec3f> positions;
    vtkm::cont::ArrayCopy(outputDataset.GetCoordinateSystem().GetData(), positions);

    // Create vector of positions to sample Velocity field and RST field
    vtkm::Id numPts = positions.GetNumberOfValues();

    // Probe velocity field at the current position
    probeFilter.SetGeometry(outputDataset);
    probeFilter.SetInvalidValue(0.0f);
    vtkm::cont::DataSet probedData = probeFilter.Execute(ds0);

    // std::cout << "Probed data size: " << probedData.GetNumberOfCells() << std::endl;

    // Get the field data from probedData at every position
    // std::cout << "Getting velocity field data... " << std::endl;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> velocityField;
    probedData.GetField("U").GetData().AsArrayHandle(velocityField);

    // Random normal vector using ArrayHandleRandomUniformReal
    vtkm::Id NumPoints = probedData.GetNumberOfCells();

    auto Normal3 = vtkm::cont::make_ArrayHandleCompositeVector(
        vtkm::cont::ArrayHandleRandomUniformReal<vtkm::FloatDefault>(NumPoints),
        vtkm::cont::ArrayHandleRandomUniformReal<vtkm::FloatDefault>(NumPoints),
        vtkm::cont::ArrayHandleRandomUniformReal<vtkm::FloatDefault>(NumPoints));

    // Random position at entrance
    std::vector<vtkm::Vec3f> randomPositionsAtEntrance_;
    for (vtkm::Id k = 0; k < NumPoints; k++) {
      randomPositionsAtEntrance_.push_back(coordinatesEntrance[rand() % coordinatesEntrance.size()]);
    }

    vtkm::cont::ArrayHandle<vtkm::Vec3f> randomPositionsAtEntrance = vtkm::cont::make_ArrayHandle(randomPositionsAtEntrance_, vtkm::CopyFlag::On);

    // std::cout << "Getting RST field data... " << std::endl;
    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 6>> RST;
    probedData.GetField("RST").GetData().AsArrayHandle(RST);

    vtkm::cont::Invoker invoke;
    auto LangevinStep = LangevinWorklet(stepSize, minCoords, maxCoords);
    invoke(LangevinStep, velocityField, RST, Normal3, randomPositionsAtEntrance, positions, flags);

    vtkm::cont::DataSetBuilderExplicitIterative dataSetBuilder;
    auto updatedPositionPortal = positions.ReadPortal();
    // NEW ----------->
    auto flagsReadPortal  = flags.ReadPortal();
    auto flagsWritePortal = flags.WritePortal();
    // <---------- NEW 
    std::vector<std::vector<double>> outputVector;

    for (vtkm::Id k = 0; k < numPts; k++) {
      vtkm::Vec3f updatedPosition = updatedPositionPortal.Get(k);
      dataSetBuilder.AddPoint(updatedPosition);
      dataSetBuilder.AddCell(vtkm::CELL_SHAPE_VERTEX);
      dataSetBuilder.AddCellPoint(k);
      outputVector.push_back({updatedPosition[0], updatedPosition[1], updatedPosition[2]});

      // NEW ----------->
      if (i % numStepsToWrite == 0) {
          rowX[k] = updatedPosition[0] - xVector[k];
          rowY[k] = updatedPosition[1] - yVector[k];
          rowZ[k] = updatedPosition[2] - zVector[k];
          rowReset[k] = flagsReadPortal.Get(k);
          flagsWritePortal.Set(k, 0); // Reset the flag to 0 after using it
      }
      // <---------- NEW 
    }

    // NEW -----------> 
    if (i % numStepsToWrite == 0) {
      dx.push_back(rowX);
      dy.push_back(rowY);
      dz.push_back(rowZ);
      spin_reset.push_back(rowReset);
    }
    // <---------- NEW

    outputDataset = dataSetBuilder.Create();
    outputVector.clear();

    std::cout << "Simulated timestep " << i << "/" << numSteps << std::endl;
  }

  Timer.Stop();
  std::cout << "Total simulation time: " << Timer.GetElapsedTime() << " seconds" << std::endl;

  // NEW -----------> Para exportar a .phantom
  Timer.Start();

  HighFive::File phantom(outputFolder + "/" + phantomFileName, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
  phantom.createAttribute<std::string>("Version", "0.9.0"); // 13/01/2025
  phantom.createAttribute<std::string>("Name", "sUbend");
  phantom.createAttribute("Ns", coordinates.size());
  phantom.createAttribute("Dims", 3);
  auto contrast = phantom.createGroup("contrast");
  contrast.createDataSet("T1", t1Vector);
  contrast.createDataSet("T2", t2Vector);
  contrast.createDataSet("ρ",  pdVector);
  auto position = phantom.createGroup("position");
  // subtractMean(xVector);
  // subtractMean(yVector);
  // subtractMean(zVector);
  position.createDataSet("x", xVector);
  position.createDataSet("y", yVector);
  position.createDataSet("z", zVector);
  auto motion   = phantom.createGroup("motion");
  auto motion_1 = motion.createGroup("motion_1");
  std::string type;
  std::string periodic;
  // Action
  auto action = motion_1.createGroup("action");
  type = "FlowPath";
  action.createAttribute<std::string>("type", type);
  action.createDataSet("dx", dx);
  action.createDataSet("dy", dy);
  action.createDataSet("dz", dz);
  action.createDataSet("spin_reset", spin_reset);
  // TimeCurve
  auto time = motion_1.createGroup("time");
  type = "TimeCurve";
  periodic = "false";
  time.createAttribute<std::string>("type", type);
  time.createAttribute<std::string>("periodic", periodic);
  time.createAttribute("periods", 1.0f);
  time.createDataSet("t", t);
  time.createDataSet("t_unit", t_unit);
  // Spins
  auto spins = motion_1.createGroup("spins");
  type = "AllSpins";
  spins.createAttribute<std::string>("type", type);

  Timer.Stop();
  std::cout << "Time spent writing .phantom file: " << Timer.GetElapsedTime() << " seconds" << std::endl;
  Timer.Reset();
  // <---------- NEW

  return 0;
}
