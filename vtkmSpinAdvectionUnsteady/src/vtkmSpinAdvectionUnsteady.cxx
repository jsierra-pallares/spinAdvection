// -----------------------------------------------------------------------------
// ██╗   ██╗████████╗██╗  ██╗███╗   ███╗        ██╗      ██████╗███████╗
// ██║   ██║╚══██╔══╝██║ ██╔╝████╗ ████║        ██║     ██╔════╝██╔════╝
// ██║   ██║   ██║   █████╔╝ ██╔████╔██║   ███  ██║     ██║     ███████╗
// ╚██╗ ██╔╝   ██║   ██╔═██╗ ██║╚██╔╝██║        ██║     ██║     ╚════██║
//  ╚████╔╝    ██║   ██║  ██╗██║ ╚═╝ ██║        ███████╗╚██████╗███████║
//   ╚═══╝     ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝        ╚══════╝ ╚═════╝╚══════╝
// -----------------------------------------------------------------------------
// Versión 0.1 - 2024, (c) José Sierra Pallares, PhD
// -----------------------------------------------------------------------------

#include <cstdlib>
#include <fstream>
#include <highfive/highfive.hpp>
#include <iostream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <vtkm/Particle.h>
#include <vtkm/Types.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/filter/flow/PathParticle.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

bool isInsideBox(const vtkm::Vec3f &position, const vtkm::Vec3f &minCoords,
                 const vtkm::Vec3f &maxCoords) {
  return (position[0] >= minCoords[0] && position[0] <= maxCoords[0] &&
          position[1] >= minCoords[1] && position[1] <= maxCoords[1] &&
          position[2] >= minCoords[2] && position[2] <= maxCoords[2]);
}

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
  auto opts = vtkm::cont::InitializeOptions::DefaultAnyDevice;
  auto config = vtkm::cont::Initialize(argc, argv, opts);

  if (argc < 1) {
    std::cerr << "Usage:\n"
              << "vtkmSpinAdvectionUnsteady settings_file.json\n"
              << std::endl;
    std::cerr << "Example:\n";
    std::cerr << "./vtkmSpinAdvectionUnsteady advection_settings.json\n";
    std::cerr << "where vtk-m options are: " << std::endl
              << config.Usage << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string settings_file = std::string(argv[1]);

  std::ifstream settingsFile(settings_file);
  if (!settingsFile.is_open()) {
    std::cerr << "Failed to open settings file." << std::endl;
    return 1;
  }

  nlohmann::json settings;
  settingsFile >> settings;

  std::string list_of_files = settings["list_of_files"];
  std::string fieldName = settings["fieldName"];
  int numStepsToWrite = settings["numStepsToWrite"];
  int numSteps = settings["numSteps"];
  std::string outputFolder = settings["outputFolder"];
  std::string seedName = settings["seedName"];
  std::string entranceName = settings["entranceName"];
  std::string solverType = settings["solverType"];
  std::string h5filename = settings["h5filename"];

  // NEW ----------->
  std::string phantomFileName = settings["phantomFileName"];
  float T1 = settings["T1"];
  float T2 = settings["T2"];
  float PD = settings["PD"];
  // <---------- NEW 

  std::cout << "fieldName: " << fieldName << std::endl;
  std::cout << "numSteps: " << numSteps << std::endl;
  std::cout << "outputFolder: " << outputFolder << std::endl;
  std::cout << "seedName: " << seedName << std::endl;
  std::cout << "entranceName: " << entranceName << std::endl;
  std::cout << "solverType: " << solverType << std::endl;
  std::cout << "h5filename: " << h5filename << std::endl;

  // NEW ----------->
  std::cout << "phantomFileName: " << phantomFileName << std::endl;
  std::cout << "T1: " << T1 << std::endl;
  std::cout << "T2: " << T2 << std::endl;
  std::cout << "PD: " << PD << std::endl;

  std::string zenodoBaseURL = "https://zenodo.org/records/16995252/files/";
  // <---------- NEW 

  // HighFive::File file(outputFolder + "/" + h5filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

  // Read coordinates from a file and store them in a vtkm::cont::ArrayHandle
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
  vtkm::Vec3f minCoords = vtkm::Vec3f(std::numeric_limits<vtkm::Float32>::max());
  vtkm::Vec3f maxCoords = vtkm::Vec3f(std::numeric_limits<vtkm::Float32>::lowest());
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

  vtkm::cont::ArrayHandle<vtkm::Particle> seeds;
  seeds.Allocate(coordinates.size());

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

  std::vector<float> t      = {0.0, 1.0};
  std::vector<float> t_unit = {0.0, 1.0};
  // <---------- NEW 

  // Write seeds to a dataset as output_0
  std::vector<std::vector<double>> outputVector;
  vtkm::cont::DataSetBuilderExplicitIterative dataSetBuilder;
  auto seedPortal = seeds.WritePortal();
  for (std::vector<vtkm::Vec3f>::size_type i = 0; i < coordinates.size(); i++) {
    vtkm::Particle p;
    p.SetPosition(coordinates[i]);
    p.SetID(i);
    seedPortal.Set(i, p);
    dataSetBuilder.AddPoint(coordinates[i]);
    dataSetBuilder.AddCell(vtkm::CELL_SHAPE_VERTEX);
    dataSetBuilder.AddCellPoint(i);
    outputVector.push_back({coordinates[i][0], coordinates[i][1], coordinates[i][2]});
    
    // NEW ----------->    
    xVector.push_back(coordinates[i][0]);
    yVector.push_back(coordinates[i][1]);
    zVector.push_back(coordinates[i][2]);
    t1Vector.push_back(T1);
    t2Vector.push_back(T2);
    pdVector.push_back(PD);
    // <---------- NEW 
  }

  // vtkm::cont::DataSet initialPositions = dataSetBuilder.Create();
  // std::string output0 = outputFolder + "/output_0.vtk";
  // vtkm::io::VTKDataSetWriter writer0(output0);
  // writer0.WriteDataSet(initialPositions);

  std::ifstream entranceSeedsFile(entranceName);
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

  // Read list of files and timestamps from file
  std::ifstream fileList(list_of_files);
  if (!fileList.is_open()) {
      std::cerr << "Failed to open file list." << std::endl;
      return 1;
  }

  std::vector<std::string> fileNames;
  std::vector<double> timestamps;
  std::string fileEntry;

  while (std::getline(fileList, fileEntry)) {
      std::istringstream iss(fileEntry);
      std::string fileName;
      double timestamp;

      if (!(iss >> fileName >> timestamp)) {
          std::cerr << "Failed to read file entry from file." << std::endl;
          return 1;
      }

      fileNames.push_back(fileName);
      timestamps.push_back(timestamp);
  }

  vtkm::filter::flow::PathParticle filter;

  if (solverType == "RK4") {
    filter.SetSolverRK4();
  } else if (solverType == "Euler") {
    filter.SetSolverEuler();
  } else {
    std::cerr << "Invalid solver type." << std::endl;
    return 1;
  }

  filter.SetUseThreadedAlgorithm(true);
  filter.SetActiveField(fieldName);
  filter.SetSeeds(seeds);

  vtkm::cont::Timer timer;

  vtkm::Float64 totalSimulationTime = 0.0;
  vtkm::Float64 totalIOTime = 0.0;

  // Store the datasets in a matrix
  timer.Start();
  std::cout << "Storing datasets in memory..." << std::endl;
  std::vector<vtkm::cont::DataSet> datasetsVector;

  // NEW ----------->
  std::vector<float> rowX(Nspins, 0.0f);
  std::vector<float> rowY(Nspins, 0.0f);
  std::vector<float> rowZ(Nspins, 0.0f);
  std::vector<int> rowReset(Nspins, 0);
  dx.push_back(rowX);
  dy.push_back(rowY);
  dz.push_back(rowZ);
  spin_reset.push_back(rowReset);

  std::vector<int> flags(Nspins, 0); // flags for every particle to store the previous "inside" boolean. Initialize to an array of zeros
  // <---------- NEW 

  #pragma omp parallel for   // Para entrada/salida con openMP. Independiente de cómo se resuelvan las trayectorias. Esto es para paralelización híbrida
  for (size_t i = 0; i < fileNames.size(); i++) {
    vtkm::io::VTKDataSetReader reader(fileNames[i]);
    vtkm::cont::DataSet dataset = reader.ReadDataSet();
    datasetsVector.push_back(dataset);
    std::cout << "Dataset " << fileNames[i] << " stored." << std::endl;
  }
  std::cout << "Datasets stored in memory." << std::endl;
  timer.Stop();

  std::cout << "Time spent storing datasets in memory: "
            << timer.GetElapsedTime() << " seconds" << std::endl;
  totalIOTime += timer.GetElapsedTime();
  timer.Reset();

  // Perform particle advection for each file and timestamp
  for (size_t i = 0; i < fileNames.size() - 1; i++) {
    timer.Start();
    vtkm::FloatDefault delta_t = std::abs(timestamps[i + 1] - timestamps[i]);
    std::cout << "Iteration number " << i << std::endl;

    vtkm::cont::DataSet ds1 = datasetsVector[i];
    vtkm::cont::DataSet ds2 = datasetsVector[i + 1];

    std::cout << "DS1: " << fileNames[i]     << ", Timestamp: " << timestamps[i]     << std::endl;
    std::cout << "DS2: " << fileNames[i + 1] << ", Timestamp: " << timestamps[i + 1] << std::endl;

    filter.SetPreviousTime(0.0);
    filter.SetNextTime(delta_t);
    filter.SetNextDataSet(ds2);
    filter.SetNumberOfSteps(numSteps);
    vtkm::Float32 stepSize = delta_t / numSteps;
    filter.SetStepSize(stepSize);

    vtkm::cont::DataSet output = filter.Execute(ds1);
    timer.Stop();

    std::cout << "Time spent per spin advection computation: " << timer.GetElapsedTime() << " seconds" << std::endl;
    totalSimulationTime += timer.GetElapsedTime();
    timer.Reset();

    // Write to VTK files
    // timer.Start();
    // std::string outputTime = outputFolder + "/output_" + std::to_string(i + 1);
    // vtkm::io::VTKDataSetWriter writer(outputTime + ".vtk");
    // writer.WriteDataSet(output);
    // timer.Stop();
    // std::cout << "Time spent writing output to disk: " << timer.GetElapsedTime() << " seconds" << std::endl;

    // Use output positions as new seed
    vtkm::cont::ArrayHandle<vtkm::Vec3f> pts;
    vtkm::cont::ArrayCopy(output.GetCoordinateSystem().GetData(), pts);
    vtkm::cont::ArrayHandle<vtkm::Particle> newSeeds;

    vtkm::Id numPts = pts.GetNumberOfValues();
    std::cout << "Number of particles: " << numPts << std::endl;
    newSeeds.Allocate(numPts);
    auto ptsPortal = pts.ReadPortal();
    auto newSeedsPortal = newSeeds.WritePortal();

    // Loop over particles to see if they have left the domain
    timer.Start();
    for (vtkm::Id j = 0; j < numPts; j++) {
      vtkm::Particle p;
      // check if positions are inside FoV
      vtkm::Vec3f position = ptsPortal.Get(j);
      // push to outputVector
      outputVector.push_back({position[0], position[1], position[2]});

      // NEW ----------->
      if (i % numStepsToWrite == 0) {
          rowX[j] = position[0] - xVector[j];
          rowY[j] = position[1] - yVector[j];
          rowZ[j] = position[2] - zVector[j];
          rowReset[j] = flags[j];
          flags[j] = 0;
      }
      // <---------- NEW 

      bool inside = isInsideBox(position, minCoords, maxCoords);
      if (!inside) {
        //  out of bounds;
        // random assignation to entrance
        vtkm::Vec3f newPosition = coordinatesEntrance[rand() % coordinatesEntrance.size()];
        p.SetID(j);
        p.SetPosition(newPosition);
        newSeedsPortal.Set(j, p);
        flags[j] = 1;
      } else {
        // still in between bounds
        p.SetID(j);
        p.SetPosition(position);
        newSeedsPortal.Set(j, p);
      }
    }

    // NEW -----------> 
    if (i % numStepsToWrite == 0) {
      dx.push_back(rowX);
      dy.push_back(rowY);
      dz.push_back(rowZ);
      spin_reset.push_back(rowReset);
    }
    // <---------- NEW

    timer.Stop();
    std::cout << "Time needed for reseeding " << timer.GetElapsedTime() << std::endl;
    totalIOTime += timer.GetElapsedTime();
    timer.Reset();

    seeds = newSeeds;
    filter.SetSeeds(seeds);

    outputVector.clear();
  }

  // NEW -----------> Para exportar a .phantom
  timer.Start();

  HighFive::File phantom(outputFolder + "/" + phantomFileName, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
  phantom.createAttribute<std::string>("Version", "0.9.0"); // 13/01/2025
  phantom.createAttribute<std::string>("Name", "Aorta");
  phantom.createAttribute("Ns", coordinates.size());
  phantom.createAttribute("Dims", 3);
  auto contrast = phantom.createGroup("contrast");
  contrast.createDataSet("T1", t1Vector);
  contrast.createDataSet("T2", t2Vector);
  contrast.createDataSet("ρ",  pdVector);
  auto position = phantom.createGroup("position");
  subtractMean(xVector);
  subtractMean(yVector);
  subtractMean(zVector);
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

  timer.Stop();
  std::cout << "Time spent writing .phantom file: " << timer.GetElapsedTime() << " seconds" << std::endl;
  timer.Reset();
  // <---------- NEW

  std::cout << "Total simulation time: " << totalSimulationTime << " seconds"
            << std::endl;

  return 0;
}
