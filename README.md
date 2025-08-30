# Installation Guide

This guide provides instructions on how to install this solver, including its dependencies:
- [VTK-m](https://gitlab.kitware.com/vtk/vtk-m) (We will upgrade soon to [Viskores](https://github.com/Viskores/viskores))
- [vcpkg](https://github.com/microsoft/vcpkg), which is included in this repository as a Git submodule. 
- HighFive, JSON, and CURL
- CMake and a C++ compatible compiler are also needed

---

## 1. Installing VTK-m

VTK-m is a toolkit for scientific visualization. Follow the instructions on the official website to install it:

- Visit the [VTK-m Readme](https://gitlab.kitware.com/vtk/vtk-m/blob/master/README.md).
- Follow the platform-specific instructions provided there.
- Be sure to compile using C++14 for compatibility. This can be done adding the appropriate flag to your compiler
- To use CUDA you need to install the CUDA toolkit from NVIDIA. Follow the instructions on the [NVIDIA CUDA Installation Guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html) for your specific operating system.
- VTK-m can be compiled for use with different accelrators. CUDA and OpenMP are supported in this project.

---
**NOTE**
Steps 2, 3, and 4 can be omitted by using:
```sh
chmod +x build.sh
./build.sh
```
---


## 2. Installing vcpkg

vcpkg is a C++ library manager that simplifies the installation of libraries. To install vcpkg:
****
1. Clone the vcpkg repository:
   ```sh
   git clone https://github.com/microsoft/vcpkg.git
   cd vcpkg
   ```

2. Bootstrap vcpkg:
   ```sh
   ./bootstrap-vcpkg.sh
   ```

3. Add vcpkg to your system's PATH (optional but recommended):
   ```sh
   ./vcpkg integrate install
   ```

For more details, visit the [vcpkg GitHub repository](https://github.com/microsoft/vcpkg).

---

## 3. Installing HDF5 and JSON Headers Using vcpkg

Once vcpkg is installed, you can use it to install the HDF5 and JSON libraries:

1. Install HighFive - HDF5 headers:
```sh
./vcpkg install highfive
```

2. Install JSON headers (e.g., nlohmann/json):
```sh
./vcpkg install nlohmann-json
```

3. Install CURL:
```sh
./vcpkg install curl
```

---
**NOTE**
If you're using a specific compiler or triplet, specify it during installation:
```sh
./vcpkg install hdf5 nlohmann-json --triplet=x64-windows
```
---


---

By following these steps, you will have VTK-m, vcpkg, and the required libraries installed and ready for use.

## 4. Compilation of the project

Source code for both solvers are provided in the repo along with a CMakeLists.txt file to facilitate compilation. Make sure you have CMake and a compatible compiler installed on your system. After successful compilation, the executables will be available in the build directory.
An example of compilation using CMake and vtk-m is provider here: [VTK-m Building Guide](https://gitlab.kitware.com/vtk/vtk-m/blob/master/README.md#building).

## 5. Minimal working examples


The source of this project includes two minimal working examples of the solvers: sUbend_mwe and aorta_mwe. Each example has a python notebbook to prepare the simulation data and scripts to run each solver.
