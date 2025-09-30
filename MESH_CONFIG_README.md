# Mesh Configuration System

## Overview

The mesh configuration is now centralized in `mesh_config.json`. This ensures that the C++ solver and Python post-processing scripts always use the same mesh parameters.

## Configuration File: `mesh_config.json`

```json
{
  "mulFac": 5,
  "meshes": {
    "bg": {
      "name": "background",
      "x0": 0.0,
      "y0": 0.0,
      "theta_deg": 0.0,
      "length": 1.0,
      "width": 1.0,
      "Nx_base": 10,
      "Ny_base": 10
    },
    "comp": {
      "name": "component",
      "x0": 0.42,
      "y0": 0.42,
      "theta_deg": 45.0,
      "length": 0.4,
      "width": 0.4,
      "Nx_base": 6,
      "Ny_base": 6
    }
  },
  "solver": {
    "Re": 100,
    "variables": ["U", "V", "P"]
  }
}
```

### Parameters

- **mulFac**: Multiplication factor for mesh resolution
- **x0, y0**: Origin coordinates of the mesh
- **theta_deg**: Rotation angle in degrees
- **length, width**: Dimensions of the mesh domain
- **Nx_base, Ny_base**: Base grid resolution (actual = base × mulFac)
- **Re**: Reynolds number (controls flow viscosity)

## Usage

### C++ Solver

The solver automatically reads `mesh_config.json`:

```cpp
#include "util/mesh_config.h"

MeshConfig config("mesh_config.json");
DataStructure bgMesh(config.getBgX0(), config.getBgY0(), config.getBgTheta(),
                     config.getBgLength(), config.getBgWidth(), 
                     config.getBgNx(), config.getBgNy());
```

### Python Post-Processing

All plotting scripts use the configuration loader:

```python
from mesh_config_loader import load_mesh_config

config = load_mesh_config()
mesh_desc = config.get_all_mesh_descriptors()
mesh_names = config.get_mesh_names()
```

## Benefits

✅ **Single Source of Truth**: Edit mesh parameters in one place
✅ **Consistency**: C++ and Python always use the same configuration
✅ **No Sync Issues**: Never worry about forgetting to update both files
✅ **Easy Experimentation**: Change parameters and rerun everything

## Workflow

1. **Edit mesh configuration**: Modify `mesh_config.json`
2. **Rebuild solver**: `cmake --build build`
3. **Run solver**: `.\build\overset.exe`
4. **Generate plots**: 
   - `python post_processing/plot_contours.py`
   - `python post_processing/plot_hole_cutting.py`
   - `python post_processing/plot_centerline.py`

All scripts automatically pick up the new configuration!

## Files

- `mesh_config.json` - Central configuration file
- `src/util/mesh_config.h` - C++ configuration reader
- `post_processing/mesh_config_loader.py` - Python configuration reader
