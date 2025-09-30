"""
Plot centerline velocity profiles and compare with Ghia et al. benchmark data.
- Reads mesh configuration from mesh_config.json
- Reads solver output from TEMP folder
- Compares U velocity at x=0.5 and V velocity at y=0.5 with Ghia et al.
- Saves PNG plots to TEMP folder

Usage (PowerShell):
  py post_processing/plot_centerline.py
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from pathlib import Path
from mesh_config_loader import load_mesh_config

DTYPE = np.float64
PI = np.pi

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 14})

# Workspace root
ROOT = Path(__file__).resolve().parents[1]

# Load mesh configuration from JSON
config = load_mesh_config()
mesh = tuple(config.get_all_mesh_descriptors())
meshNames = config.get_mesh_names()
varNames = config.get_var_names()
Re = config.solver.get('Re', 100)  # Get Reynolds number from config

# Ghia et al. benchmark data (Re=100, 400, and 1000)
GHIA_DATA = {
    100: {
        "U": {
            "y": [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 
                  0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000],
            "u": [1.00000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, 
                  -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.00000]
        },
        "V": {
            "x": [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 
                  0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000],
            "v": [0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 
                  0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000]
        }
    },
    400: {
        "U": {
            "y": [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 
                  0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000],
            "u": [1.00000, 0.75837, 0.68439, 0.61756, 0.55892, 0.29093, 0.16256, 0.02135, 
                  -0.11477, -0.17119, -0.32726, -0.24299, -0.14612, -0.10338, -0.09266, -0.08186, 0.00000]
        },
        "V": {
            "x": [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 
                  0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000],
            "v": [0.00000, -0.12146, -0.15663, -0.19254, -0.22847, -0.23827, -0.44993, -0.38598, 
                  0.05186, 0.30174, 0.30203, 0.28124, 0.22965, 0.20920, 0.19713, 0.18360, 0.00000]
        }
    },
    1000: {
        "U": {
            "y": [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 
                  0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000],
            "u": [1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, 
                  -0.06080, -0.10648, -0.27805, -0.38289, -0.29730, -0.22220, -0.20196, -0.18109, 0.00000]
        },
        "V": {
            "x": [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 
                  0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000],
            "v": [0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.51500, -0.42665, -0.31966, 
                  0.02526, 0.32235, 0.33075, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0.00000]
        }
    }
}


def plot_ghia(varName: str, Re: int):
    """Plot Ghia et al. benchmark data for the specified Reynolds number."""
    # Find closest available Re in Ghia data
    available_Re = sorted(GHIA_DATA.keys())
    closest_Re = min(available_Re, key=lambda x: abs(x - Re))
    
    if abs(closest_Re - Re) > 50:
        print(f"  Warning: No exact Ghia data for Re={Re}, using closest Re={closest_Re}")
    
    ghia_data = GHIA_DATA[closest_Re]
    
    if varName == "U":
        axis = ghia_data["U"]["y"]
        var = ghia_data["U"]["u"]
    elif varName == "V":
        axis = ghia_data["V"]["x"]
        var = ghia_data["V"]["v"]
    else:
        return
    
    plt.scatter(axis, var, c="black", marker="^", s=50, 
               label=f"Ghia et al. (Re={closest_Re})", zorder=10)


def centerline(title: str, filename: str, velComp: str, output_dir: Path):
    """Generate centerline plots for specified velocity component."""
    print(f"Plotting - {velComp}")
    
    fig, ax = plt.subplots(figsize=(9, 8))
    centerN = 100

    # Define centerline coordinates
    if velComp == "U":
        xCenterAxis = np.ones(1) * 0.5
        yCenterAxis = np.linspace(0, 1, centerN, endpoint=True)
    elif velComp == "V":
        xCenterAxis = np.linspace(0, 1, centerN, endpoint=True)
        yCenterAxis = np.ones(1) * 0.5
    else:
        print(f"Warning: Unknown velocity component '{velComp}'")
        return

    xCenterGrid, yCenterGrid = np.meshgrid(xCenterAxis, yCenterAxis)
    k = varNames.index(velComp)
    
    # Plot Ghia benchmark first
    plot_ghia(velComp, Re)

    # Process each mesh
    for idx, mesh_name in enumerate(meshNames):
        meshDetails = mesh[idx]
        xOrigin, yOrigin, Theta, L, W, Nx, Ny = list(meshDetails)
        Nx, Ny = int(Nx) + 2, int(Ny) + 2
        print(f"  Mesh '{mesh_name}' discretization: {Nx} x {Ny}, {Theta * 180.0/PI:.1f} deg")

        # Load solver output data
        dataFile = output_dir / f"output_Upwind_{mesh_name}Mesh.dat"
        if not dataFile.exists():
            print(f"  Warning: {dataFile} not found. Skipping {mesh_name} mesh.")
            continue
            
        varData = np.loadtxt(dataFile, comments='#', dtype=DTYPE)
        varData = varData[k*Ny:(k+1)*Ny, :]

        # Generate mesh grid
        xAxis = np.linspace(0, L, Nx) 
        yAxis = np.linspace(0, W, Ny)
        _xGrid, _yGrid = np.meshgrid(xAxis, yAxis)
        _xGrid = _xGrid + L/Nx/2
        _yGrid = _yGrid + W/Ny/2
        xGrid = xOrigin + _xGrid * np.cos(Theta) - _yGrid * np.sin(Theta)
        yGrid = yOrigin + _xGrid * np.sin(Theta) + _yGrid * np.cos(Theta)

        # Load point type mask
        maskFile = output_dir / f"{mesh_name}PtType.txt"
        if not maskFile.exists():
            print(f"  Warning: {maskFile} not found. Skipping {mesh_name} mesh.")
            continue
            
        pointMask = np.loadtxt(maskFile)
        cellTypeMask = np.where((pointMask == 0), False, True)
        field_var = varData.flatten()
        field_var[np.invert(cellTypeMask)] = np.nan

        # Interpolate to centerline
        interpolatedCenterLine = griddata(
            np.stack((xGrid.flatten(), yGrid.flatten()), axis=1), 
            field_var, 
            (xCenterGrid, yCenterGrid), 
            method='linear'
        )

        # Plot centerline
        if velComp == "U":
            ax.plot(yCenterGrid.flatten(), interpolatedCenterLine.flatten(), 
                   label=f"{mesh_name} Mesh", linewidth=2)
        elif velComp == "V":
            ax.plot(xCenterGrid.flatten(), interpolatedCenterLine.flatten(), 
                   label=f"{mesh_name} Mesh", linewidth=2)

    # Finalize plot
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    fig.suptitle(f"{title} for {velComp}")
    ax.set_ylabel(velComp)
    
    if velComp == "U":
        ax.set_xlabel("y")
    elif velComp == "V":
        ax.set_xlabel("x")
    
    # Save figure
    out_path = output_dir / f"{filename}_{velComp}.png"
    fig.tight_layout()
    plt.savefig(out_path, bbox_inches="tight", pad_inches=0.2, dpi=150)
    plt.close()
    print(f"Saved: {out_path}")


def main(output_dir: Path = None):
    if output_dir is None:
        # Default to latest directory in TEMP
        temp_dir = ROOT / "TEMP"
        output_dirs = [d for d in temp_dir.iterdir() if d.is_dir() and d.name.startswith("Re")]
        if not output_dirs:
            raise FileNotFoundError("No output directories found. Run solver first!")
        output_dir = max(output_dirs, key=lambda d: d.stat().st_mtime)
    
    output_dir = Path(output_dir)
    print("Generating centerline plots...")
    title = "Centerline plot"
    filename = "lid_overset_centerline"
    
    for var in varNames[:2]:  # U and V only
        centerline(title, filename, var, output_dir)
    
    print("Done!")


if __name__ == "__main__":
    main()
