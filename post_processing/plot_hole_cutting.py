"""
Visualize hole cutting results showing point types for bg/comp meshes.
- Reads mesh configuration from mesh_config.json
- Reads point type files from TEMP folder
- Saves PNG plots to TEMP folder

Point types:
  0: unused
  1: calculated
  2: donor
  3: receiver
  4: donor buffer
  5: bc specified

Usage (PowerShell):
  py post_processing/plot_hole_cutting.py
"""
from __future__ import annotations
import matplotlib as mpl
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
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
mesh = np.array(config.get_all_mesh_descriptors())
meshNames = config.get_mesh_names()

# Point type configuration
pointTypeLabels = ["unused", "calculated", "donor", "receiver", "donor buffer", "bc specified"]
pointType = np.linspace(-0.5, 5.5, len(pointTypeLabels)+1)
pointTypeLabels = [str(int(m))+" : "+n for m, n in zip(pointType+0.5, pointTypeLabels)]

def plot_hole_cutting(output_dir: Path):
    """Generate hole cutting visualization for both meshes."""
    cmap = mpl.colormaps["turbo"]
    norm = colors.BoundaryNorm(pointType, cmap.N)

    for meshIdx in range(len(meshNames)):
        meshDetails = mesh[meshIdx] 
        xOrigin, yOrigin = meshDetails[0], meshDetails[1]
        theta = meshDetails[2]
        L, W = meshDetails[3], meshDetails[4]
        Nx, Ny = int(meshDetails[5]), int(meshDetails[6])
        
        # Generate mesh grid (with halo: +2)
        _x = np.linspace(0, L, Nx+2)
        _y = np.linspace(0, W, Ny+2)
        xLocal, yLocal = np.meshgrid(_x, _y)
        
        # Rotate and translate
        x = xOrigin + xLocal*np.cos(theta) - yLocal*np.sin(theta)
        y = yOrigin + xLocal*np.sin(theta) + yLocal*np.cos(theta)
        
        # Load point type data
        ptTypeFile = output_dir / f"{meshNames[meshIdx]}PtType.txt"
        if not ptTypeFile.exists():
            print(f"Warning: {ptTypeFile} not found. Skipping {meshNames[meshIdx]} mesh.")
            continue
            
        pointTypeData = np.loadtxt(ptTypeFile)
        print(f"Mesh '{meshNames[meshIdx]}': Point type data shape = {pointTypeData.shape}")
        
        # Calculate image scale based on component mesh
        imageScale = np.floor([Nx/mesh[1][-2], Ny/mesh[1][-1]])
        
        # Create figure
        fig, ax = plt.subplots(figsize=(7 * imageScale[0], 5 * imageScale[1]))
        im = ax.scatter(x, y, c=pointTypeData, cmap=cmap, norm=norm, s=20)
        
        # Add colorbar
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_ticks(ticks=pointType[:-1]+0.5, labels=pointTypeLabels)
        
        # Labels
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(f'Hole Cutting - {meshNames[meshIdx].upper()} Mesh')
        ax.set_aspect('equal')
        
        # Save figure
        out_path = output_dir / f"{meshNames[meshIdx]}HoleCut.png"
        fig.savefig(out_path, bbox_inches="tight", pad_inches=0.2, dpi=150)
        plt.close(fig)
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
    print("Generating hole cutting visualizations...")
    plot_hole_cutting(output_dir)
    print("Done!")


if __name__ == "__main__":
    main()
