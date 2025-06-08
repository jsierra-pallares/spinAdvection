rm -r ./output/*

vtkmSpinAdvectionUnsteady  --vtkm-device="OpenMP" advection_settings.json 
