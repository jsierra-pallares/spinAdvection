#!/bin/bash

# Initialize vcpkg
cd vcpkg
./bootstrap-vcpkg.sh
cd ..

# vtkmSpinAdvectionUnsteady
echo "------------- Building vtkmSpinAdvectionUnsteady -------------"
cd vtkmSpinAdvectionUnsteady
mkdir build
cmake -B build -S .
cmake --build build
cmake --install build
export PATH=$PATH:$PWD/build
cd ..

# vtkmTurbulentTrajectory
echo "------------- Building vtkmTurbulentTrajectory -------------"
cd vtkmTurbulentTrajectory
mkdir build
cmake -B build -S .
cmake --build build
cmake --install build
export PATH=$PATH:$PWD/build
cd ..