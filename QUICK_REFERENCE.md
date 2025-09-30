# Quick Reference Guide

## Complete Workflow

### 1. Configure Parameters
Edit `mesh_config.json`:
```json
{
  "mulFac": 5,
  "meshes": {
    "bg": { "x0": 0.0, "y0": 0.0, "theta_deg": 0.0, ... },
    "comp": { "x0": 0.42, "y0": 0.42, "theta_deg": 45.0, ... }
  },
  "solver": {
    "Re": 100
  }
}
```

### 2. Build & Run Solver
```powershell
cmake --build build
.\build\overset.exe
```

### 3. Generate All Plots
```powershell
python post_processing\plot_all.py
```

### 4. View Results
All output files are in `TEMP/`:
- Solver data: `output_Upwind_*.dat`, `*PtType.txt`
- Plots: `*.png` files

---

## File Structure

```
mesh_config.json                    # Configuration file (EDIT THIS!)
build/overset.exe                   # Compiled solver
TEMP/                               # All outputs go here
  ├── output_Upwind_bgMesh.dat      # Solver output
  ├── output_Upwind_compMesh.dat
  ├── bgPtType.txt                  # Point type masks
  ├── compPtType.txt
  ├── lid_overset_U.png             # Contour plots
  ├── lid_overset_V.png
  ├── lid_overset_P.png
  ├── bgHoleCut.png                 # Hole cutting
  ├── compHoleCut.png
  ├── lid_overset_centerline_U.png  # Centerline plots
  └── lid_overset_centerline_V.png
post_processing/
  ├── plot_all.py                   # Run all plots (USE THIS!)
  ├── plot_contours.py              # Individual scripts
  ├── plot_hole_cutting.py
  ├── plot_centerline.py
  └── mesh_config_loader.py         # Config reader
```

---

## Common Tasks

### Change Reynolds Number
1. Edit `mesh_config.json`: `"Re": 1000`
2. Rebuild: `cmake --build build`
3. Run: `.\build\overset.exe`

### Change Mesh Rotation
1. Edit `mesh_config.json`: `"theta_deg": 30.0`
2. Rebuild: `cmake --build build`
3. Run: `.\build\overset.exe`

### Change Grid Resolution
1. Edit `mesh_config.json`: `"mulFac": 10`
2. Rebuild: `cmake --build build`
3. Run: `.\build\overset.exe`

### Re-generate Plots Only
```powershell
python post_processing\plot_all.py
```

---

## Troubleshooting

### Build Errors
- Check MinGW is in PATH: `$env:Path = "C:\MinGW\bin;" + $env:Path`
- Clean build: `Remove-Item -Recurse build/CMakeFiles; cmake --build build`

### Python Errors
- Install dependencies: `pip install numpy matplotlib scipy`
- Check files exist: `ls TEMP\*.dat`

### Config Not Working
- After editing `mesh_config.json`, always rebuild: `cmake --build build`
- The C++ code must be recompiled to read the new config

---

## Output Files Description

| File | Description |
|------|-------------|
| `output_Upwind_*Mesh.dat` | Flow field data (U, V, P) |
| `*PtType.txt` | Point type masks (donor, receiver, etc.) |
| `lid_overset_*.png` | Main contour visualizations |
| `*HoleCut.png` | Hole cutting visualization |
| `lid_overset_centerline_*.png` | Comparison with Ghia et al. benchmark |

---

## Tips

✅ **Single command workflow:**
```powershell
cmake --build build; .\build\overset.exe; python post_processing\plot_all.py
```

✅ **Keep TEMP clean:** Output files are isolated in TEMP folder, safe to delete anytime

✅ **Version control:** TEMP/ is in `.gitignore`, won't clutter your commits

✅ **Experimentation:** Change `mesh_config.json` parameters freely, everything syncs automatically
