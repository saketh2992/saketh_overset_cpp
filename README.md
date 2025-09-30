# Overset Mesh Solver (C++)

A simple overset mesh CFD solver written in C++ with Python post-processing tools.

## Quick Start

1. **Configure mesh parameters** in `mesh_config.json` (mesh size, rotation, Reynolds number)
2. **Build and run solver**: `cmake --build build; .\build\overset.exe`
3. **Generate all plots**: `python post_processing\plot_all.py`

All output files are saved to the `TEMP/` folder.

## Configuration

Edit `mesh_config.json` to set:
- Mesh positions, rotations, and grid resolution
- Reynolds number (Re)
- Other solver parameters

See `MESH_CONFIG_README.md` for details.

## Build (MinGW + CMake)

```powershell
# Ensure MinGW is on PATH
$env:Path = "C:\\MinGW\\bin;" + $env:Path

# Configure & build
cmake --build build

# Run
.\build\overset.exe
```

## Post-Processing

After a successful run, all output data is written to the `TEMP/` folder.

### Generate All Plots (Recommended)

Run a single command to generate all visualizations:

```powershell
# Optional: create a virtual environment and install dependencies (one-time)
py -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install numpy matplotlib scipy

# Generate all plots
py post_processing\plot_all.py
```

This creates:
- **Contour plots**: `lid_overset_U.png`, `lid_overset_V.png`, `lid_overset_P.png`
- **Hole cutting**: `bgHoleCut.png`, `compHoleCut.png`
- **Centerline plots**: `lid_overset_centerline_U.png`, `lid_overset_centerline_V.png`

### Individual Plot Scripts

You can also run individual scripts:

```powershell
py post_processing\plot_contours.py      # Contour plots
py post_processing\plot_hole_cutting.py  # Hole cutting visualization
py post_processing\plot_centerline.py    # Centerline comparison with Ghia et al.
```

### Alternative: Gnuplot (no Python required)

The solver also writes gnuplot-friendly data and script to `TEMP/plot_contours.plt`:

```powershell
# Run solver first
.\build\overset.exe
# Option A: winget (Windows 10/11)
winget install Gnuplot.Gnuplot
# Option B: Chocolatey
# choco install gnuplot

# Generate PNG contours
gnuplot plot_contours.plt
```

## Notes
- Requires MinGW (g++, gcc, mingw32-make).
- CMake presets are provided in `CMakePresets.json`.
- Source code is C++11-compatible (works with GCC 6.x).
