# rm -r ./output/*

# vtkmTurbulentTrajectory --vtkm-device="OpenMP" advection_settings.json
vtkmTurbulentTrajectory --vtkm-device="Cuda" advection_settings.json
 