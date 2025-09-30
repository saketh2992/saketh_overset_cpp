"""
Plot U, V, P contours from solver outputs for bg/comp meshes.
- Reads mesh configuration from mesh_config.json
- Expects files in TEMP folder:
  - output_Upwind_bgMesh.dat
  - output_Upwind_compMesh.dat
  - bgPtType.txt
  - compPtType.txt
- Saves PNGs into TEMP/ as lid_overset_{U,V,P}.png

Usage (PowerShell):
  # Optional: create venv & install deps
  # py -m venv .venv; .\.venv\Scripts\Activate.ps1; pip install numpy matplotlib
  py post_processing/plot_contours.py
"""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from pathlib import Path
from mesh_config_loader import load_mesh_config

DTYPE = np.float64
PI = np.pi
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({"font.size": 14})

# Workspace root assumed as this file's parent directory
ROOT = Path(__file__).resolve().parents[1]

# Load mesh configuration from JSON
config = load_mesh_config()
mesh_desc = config.get_all_mesh_descriptors()
mesh_names = config.get_mesh_names()
var_names = config.get_var_names()
colours = ["k", "r", "b"]
DEFAULT_LEVELS = 20


def load_field_block(fp: Path, k: int, Ny: int) -> np.ndarray:
    """Load block k from the ASCII table where data is arranged vertically by blocks of Ny rows.
    The files contain header lines starting with '#', which loadtxt can skip via comments.
    """
    arr = np.loadtxt(fp, comments="#", dtype=DTYPE)
    # Defensive: if arr is 1D (single row), reshape to (1, -1)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    start = k * Ny
    end = (k + 1) * Ny
    if end > arr.shape[0]:
        raise ValueError(f"File {fp.name}: requested block k={k} with Ny={Ny} exceeds rows {arr.shape[0]}")
    return arr[start:end, :]


def field_var_contour(title: str, filename_prefix: str, k: int, output_dir: Path) -> None:
    print(f"Plotting - {var_names[k]}")
    
    # Files
    data_files = {
        "bg": output_dir / "output_Upwind_bgMesh.dat",
        "comp": output_dir / "output_Upwind_compMesh.dat",
    }
    mask_files = {
        "bg": output_dir / "bgPtType.txt",
        "comp": output_dir / "compPtType.txt",
    }

    fig, ax = plt.subplots(figsize=(9, 8))
    ax.tick_params(left=False, bottom=False, labelbottom=False, labelleft=False)

    var_min = None
    var_max = None
    n_levels = DEFAULT_LEVELS

    for idx, mesh_name in enumerate(mesh_names):
        x0, y0, theta, L, W, Nx, Ny = mesh_desc[idx]
        Nx, Ny = int(Nx) + 2, int(Ny) + 2
        print(f"Mesh '{mesh_name}' discretization- {Nx} {Ny} {theta * 180.0/PI:.1f} deg")

        data_fp = data_files[mesh_name]
        mask_fp = mask_files[mesh_name]
        if not data_fp.exists():
            raise FileNotFoundError(f"Missing data file: {data_fp}")
        if not mask_fp.exists():
            raise FileNotFoundError(f"Missing mask file: {mask_fp}")

        field_block = load_field_block(data_fp, k=k, Ny=Ny)

        # cell-centered uniform grid (rotate by theta)
        x_axis = np.linspace(0.0, L, Nx)
        y_axis = np.linspace(0.0, W, Ny)
        _x, _y = np.meshgrid(x_axis, y_axis)
        _x = _x + L / Nx / 2.0
        _y = _y + W / Ny / 2.0
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        xg = x0 + _x * cos_t - _y * sin_t
        yg = y0 + _x * sin_t + _y * cos_t

        # mask invalid points
        point_mask = np.loadtxt(mask_fp, dtype=DTYPE)
        xg = xg.flatten()
        yg = yg.flatten()
        fv = field_block.flatten()
        # mask: keep where point_mask != 0
        cell_mask = np.where(point_mask.flatten() == 0, False, True)
        if cell_mask.shape[0] != xg.shape[0]:
            # If mask size mismatches, attempt to reshape based on Ny,Nx, else raise
            try:
                cell_mask = point_mask.reshape(Ny, Nx).flatten() != 0
            except Exception as e:
                raise ValueError(
                    f"Mask size mismatch for {mesh_name}: mask shape {point_mask.shape}, expected {(Ny, Nx)}"
                ) from e
        xg = xg[cell_mask]
        yg = yg[cell_mask]
        fv = fv[cell_mask]

        triang = tri.Triangulation(xg, yg)
        x_rem = xg[triang.triangles] - np.roll(xg[triang.triangles], 1, axis=1)
        y_rem = yg[triang.triangles] - np.roll(yg[triang.triangles], 1, axis=1)
        maxi = np.max(np.sqrt(x_rem**2 + y_rem**2), axis=1).reshape(-1)
        triang.set_mask(maxi > np.hypot(L / Nx, W / Ny) * 1.1)

        # compute global min/max on first mesh for consistent levels
        if idx == 0:
            sorted_unique = np.sort(np.unique(fv))
            vmax = sorted_unique[-5] if sorted_unique.size >= 5 else sorted_unique[-1]
            vmin = fv.min()
            if k == 2:  # pressure, give more contour lines
                vmax = 0.8 * vmax
                n_levels = DEFAULT_LEVELS * 5
            var_min, var_max = vmin, vmax
            print(f"value range: min={var_min:.6g}, max={var_max:.6g}, levels={n_levels}")

        levels = np.linspace(var_min, var_max, n_levels)
        cs = ax.tricontour(triang, fv, levels=levels, colors=colours[idx % len(colours)])
        ax.clabel(cs, inline=1, fontsize=10)

        # draw the mesh rectangle outline
        rect = plt.Rectangle((x0 + L / Nx / 2.0, y0 + W / Ny / 2.0), L, W,
                             angle=theta * 180.0 / PI, linewidth=2, facecolor="none",
                             edgecolor=colours[idx % len(colours)])
        ax.add_artist(rect)

    fig.suptitle(f"{title} - {var_names[k]}")
    out_path = output_dir / f"{filename_prefix}_{var_names[k]}.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
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
    title = "NS contour for 2d plate"
    fname = "lid_overset"
    for k in range(3):
        field_var_contour(title, fname, k, output_dir)


if __name__ == "__main__":
    main()
