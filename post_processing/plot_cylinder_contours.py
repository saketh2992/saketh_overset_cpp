"""
Plot U, V, P contours for cylinder flow simulation.
- Reads output files from TEMP/CylinderFlow_Re100_bg100x100_cyl50x50_R0.50/
- Plots background rectangular mesh (contains the actual flow solution)
- Adds cylinder outline for visualization
- Saves PNGs as cylinder_flow_{U,V,P}.png

Usage (PowerShell):
  py post_processing/plot_cylinder_contours.py
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Configuration
ROOT = Path(__file__).resolve().parents[1]
TEMP_DIR = ROOT / "TEMP"

# Cylinder flow parameters (from main_cylinder_flow.cpp)
BG_XMIN, BG_XMAX = -5.0, 10.0
BG_YMIN, BG_YMAX = -5.0, 5.0
BG_NX, BG_NY = 50, 50  # Updated to match current simulation
CYLINDER_CENTER = (0.0, 0.0)
CYLINDER_RADIUS = 0.5
COMP_OUTER_RADIUS = 2.0  # Component mesh outer radius

VAR_NAMES = ["U", "V", "P"]
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({"font.size": 12})


def load_field_from_bgmesh(filepath, nx, ny, nvar=3):
    """
    Load flow field from background mesh output file.
    File format: blocks of data for each variable k, with (ny+2) rows and (nx+2) columns per row
    Each row contains all nx+2 values for a given j index
    """
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Split by the header markers to separate variables
    blocks = content.split('########## Data for k =')
    
    fields = []
    ny_total = ny + 2
    nx_total = nx + 2
    
    for k in range(nvar):
        if k + 1 >= len(blocks):
            print(f"Warning: Could not find data for variable k={k}")
            break
        
        # Extract the data section for this variable
        block_text = blocks[k + 1]
        
        # Parse all numbers from this block
        numbers = []
        for line in block_text.split('\n'):
            line = line.strip()
            if line and not line.startswith('#'):
                # Split by whitespace and convert to float
                values = line.split()
                for val in values:
                    try:
                        numbers.append(float(val))
                    except ValueError:
                        pass  # Skip non-numeric tokens
        
        # Reshape into (ny+2, nx+2) matrix
        expected_size = ny_total * nx_total
        if len(numbers) >= expected_size:
            field_data = np.array(numbers[:expected_size]).reshape(ny_total, nx_total)
            fields.append(field_data)
            print(f"Loaded variable {k} ({['U', 'V', 'P'][k]}): shape={field_data.shape}, "
                  f"range=[{field_data.min():.4f}, {field_data.max():.4f}]")
        else:
            print(f"Warning: Variable k={k} has {len(numbers)} numbers, expected {expected_size}")
    
    return fields


def plot_contour(x, y, field, var_name, cylinder_center, cylinder_radius, output_path):
    """
    Plot a single contour with cylinder outline
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Create contour plot
    levels = 20
    if var_name == "P":
        levels = 50  # More levels for pressure
    
    # Plot filled contours
    cf = ax.contourf(x, y, field, levels=levels, cmap='RdBu_r')
    plt.colorbar(cf, ax=ax, label=var_name)
    
    # Add contour lines
    cs = ax.contour(x, y, field, levels=10, colors='black', linewidths=0.5, alpha=0.4)
    ax.clabel(cs, inline=True, fontsize=8, fmt='%0.2f')
    
    # Draw cylinder
    circle = plt.Circle(cylinder_center, cylinder_radius, 
                       color='white', fill=True, linewidth=2, edgecolor='black', zorder=10)
    ax.add_patch(circle)
    
    # Add cylinder label
    ax.text(cylinder_center[0], cylinder_center[1], 'Cylinder\n(No-slip)', 
           ha='center', va='center', fontsize=10, weight='bold', zorder=11)
    
    # Formatting
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Flow Past Circular Cylinder - {var_name} (Re=100)')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Save
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_path}")


def main():
    # Find cylinder flow directory
    output_dirs = [d for d in TEMP_DIR.iterdir() 
                  if d.is_dir() and d.name.startswith("CylinderFlow")]
    
    if not output_dirs:
        print("Error: No CylinderFlow output directory found in TEMP/")
        print("Please run the overset_cylinder.exe first!")
        return
    
    # Use the most recently modified directory
    output_dir = max(output_dirs, key=lambda d: d.stat().st_mtime)
    print(f"Reading from: {output_dir}")
    
    # Load background mesh data
    bg_file = output_dir / "output_bgMesh.dat"
    
    if not bg_file.exists():
        print(f"Error: {bg_file} not found!")
        return
    
    print(f"Loading background mesh data from {bg_file.name}...")
    fields = load_field_from_bgmesh(bg_file, BG_NX, BG_NY, nvar=3)
    
    if len(fields) < 3:
        print("Error: Could not load all 3 fields (U, V, P)")
        return
    
    # Create coordinate grid (cell centers)
    dx = (BG_XMAX - BG_XMIN) / BG_NX
    dy = (BG_YMAX - BG_YMIN) / BG_NY
    
    # Grid for plotting (interior cells + ghost cells)
    x = np.linspace(BG_XMIN - dx/2, BG_XMAX + dx/2, BG_NX + 2)
    y = np.linspace(BG_YMIN - dy/2, BG_YMAX + dy/2, BG_NY + 2)
    X, Y = np.meshgrid(x, y)
    
    # Load point type mask (0 = unused/hole, others = active)
    pttype_file = output_dir / "bgPtType.txt"
    if pttype_file.exists():
        pt_mask = np.loadtxt(pttype_file)
        pt_mask = pt_mask.reshape(BG_NY + 2, BG_NX + 2)
        
        # Mask out hole-cut regions
        for i, field in enumerate(fields):
            fields[i] = np.ma.masked_where(pt_mask == 0, field)
    
    print(f"Grid dimensions: {BG_NX+2} x {BG_NY+2}")
    print(f"U range: [{fields[0].min():.4f}, {fields[0].max():.4f}]")
    print(f"V range: [{fields[1].min():.4f}, {fields[1].max():.4f}]")
    print(f"P range: [{fields[2].min():.4f}, {fields[2].max():.4f}]")
    
    # Plot each variable
    print("\nGenerating contour plots...")
    for k, var_name in enumerate(VAR_NAMES):
        output_path = output_dir / f"cylinder_flow_{var_name}.png"
        plot_contour(X, Y, fields[k], var_name, CYLINDER_CENTER, CYLINDER_RADIUS, output_path)
    
    print("\n✓ All contour plots generated successfully!")
    print(f"  Location: {output_dir}")


if __name__ == "__main__":
    main()
