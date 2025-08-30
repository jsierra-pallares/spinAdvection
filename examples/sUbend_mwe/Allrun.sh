OUT_DIR="output"

if [ -d "$OUT_DIR" ]; then
    rm -r "$OUT_DIR"/*
else
    mkdir "$OUT_DIR"
fi

vtkmTurbulentTrajectory --vtkm-device="Cuda" advection_settings.json
