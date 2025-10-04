"""
Visualize hole cutting for cylinder flow simulation.
- Shows point types for background and circular component meshes
- Displays how the cylinder creates a hole in the background mesh
- Saves PNG plots showing the overset mesh connectivity

Point types:
  0: UNUSED (hole-cut region inside cylinder)
  1: CALCULATED (active interior points)
  2: DONOR (provides interpolation data)
  3: INTERPOLATION_RECIEVER (boundary receiving interpolation)
  4: DONOR_BUFFER
  5: BC_SPECIFIED

Usage (PowerShell):
  py post_processing/plot_hole_cutting_cylinder.py
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors
import numpy as np
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
COMP_OUTER_RADIUS = 2.0  # Updated to match current simulation
COMP_NTHETA = 25  # Updated to match current simulation
COMP_NR = 25  # Updated to match current simulation

# Point type labels
POINT_TYPE_LABELS = {
    0: "UNUSED (hole)",
    1: "CALCULATED",
    2: "DONOR",
    3: "INTERPOLATION_RECIEVER",
    4: "DONOR_BUFFER",
    5: "BC_SPECIFIED"
}

# Color map for point types
POINT_TYPE_COLORS = ['gray', 'lightblue', 'orange', 'red', 'yellow', 'green']

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 11})


def plot_background_hole_cutting(output_dir):
    """Plot hole cutting visualization for background rectangular mesh."""
    print("\n=== Background Mesh Hole Cutting ===")
    
    # Load point type data
    pt_file = output_dir / "bgPtType.txt"
    if not pt_file.exists():
        print(f"Error: {pt_file} not found!")
        return
    
    point_types = np.loadtxt(pt_file)
    print(f"Point type data shape: {point_types.shape}")
    
    # Reshape to grid (NY+2, NX+2) - includes ghost cells
    ny_total, nx_total = BG_NY + 2, BG_NX + 2
    if point_types.size == ny_total * nx_total:
        point_types = point_types.reshape(ny_total, nx_total)
    else:
        print(f"Warning: Expected {ny_total * nx_total} points, got {point_types.size}")
        return
    
    # Create coordinate grid
    dx = (BG_XMAX - BG_XMIN) / BG_NX
    dy = (BG_YMAX - BG_YMIN) / BG_NY
    x = np.linspace(BG_XMIN - dx/2, BG_XMAX + dx/2, nx_total)
    y = np.linspace(BG_YMIN - dy/2, BG_YMAX + dy/2, ny_total)
    X, Y = np.meshgrid(x, y)
    
    # Count point types
    print("\nPoint type distribution:")
    for pt_val, pt_label in POINT_TYPE_LABELS.items():
        count = np.sum(point_types == pt_val)
        percentage = 100.0 * count / point_types.size
        print(f"  {pt_val}: {pt_label:25s} - {count:5d} points ({percentage:5.2f}%)")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Create custom colormap
    cmap = colors.ListedColormap(POINT_TYPE_COLORS)
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    
    # Plot as scatter to see individual points
    scatter = ax.scatter(X.flatten(), Y.flatten(), c=point_types.flatten(), 
                        cmap=cmap, norm=norm, s=30, marker='s', edgecolors='none')
    
    # Add colorbar
    cbar = fig.colorbar(scatter, ax=ax, ticks=[0, 1, 2, 3, 4, 5])
    cbar.ax.set_yticklabels([POINT_TYPE_LABELS[i] for i in range(6)])
    
    # Draw cylinder outline
    cylinder = plt.Circle(CYLINDER_CENTER, CYLINDER_RADIUS, 
                         fill=False, edgecolor='black', linewidth=3, linestyle='--',
                         label=f'Cylinder (R={CYLINDER_RADIUS})')
    ax.add_patch(cylinder)
    
    # Draw component mesh outer boundary
    comp_outer = plt.Circle(CYLINDER_CENTER, COMP_OUTER_RADIUS, 
                           fill=False, edgecolor='purple', linewidth=2, linestyle=':',
                           label=f'Component Outer (R={COMP_OUTER_RADIUS})')
    ax.add_patch(comp_outer)
    
    # Formatting
    ax.set_xlabel('x (m)', fontsize=12)
    ax.set_ylabel('y (m)', fontsize=12)
    ax.set_title('Background Mesh - Hole Cutting Around Cylinder', fontsize=14, weight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    ax.legend(loc='upper right', fontsize=10)
    
    # Set limits with some margin
    margin = 1.0
    ax.set_xlim(BG_XMIN - margin, BG_XMAX + margin)
    ax.set_ylim(BG_YMIN - margin, BG_YMAX + margin)
    
    # Save figure
    out_path = output_dir / "bgHoleCut.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\n✓ Saved: {out_path}")


def plot_circular_mesh_connectivity(output_dir):
    """Plot circular component mesh showing point types."""
    print("\n=== Circular Component Mesh ===")
    
    # Load point type data
    pt_file = output_dir / "compPtType.txt"
    if not pt_file.exists():
        print(f"Error: {pt_file} not found!")
        return None, None, None
    
    point_types = np.loadtxt(pt_file)
    print(f"Point type data shape: {point_types.shape}")
    
    # Infer grid dimensions from actual data size
    total_points = point_types.size
    # Try to find factors close to expected COMP_NR and COMP_NTHETA
    possible_dims = []
    for nr in range(COMP_NR - 2, COMP_NR + 3):
        for nt in range(COMP_NTHETA - 2, COMP_NTHETA + 3):
            if nr * nt == total_points:
                possible_dims.append((nt, nr))
    
    if not possible_dims:
        # Fallback: try perfect square root
        sqrt_n = int(np.sqrt(total_points))
        if sqrt_n * sqrt_n == total_points:
            ntheta_total, nr_total = sqrt_n, sqrt_n
        else:
            print(f"Error: Cannot determine grid dimensions for {total_points} points")
            return None, None, None
    else:
        ntheta_total, nr_total = possible_dims[0]
    
    print(f"Using grid dimensions: {ntheta_total} (theta) x {nr_total} (radial)")
    point_types = point_types.reshape(ntheta_total, nr_total)
    
    # Create polar grid
    r = np.linspace(CYLINDER_RADIUS, COMP_OUTER_RADIUS, nr_total)
    theta = np.linspace(0, 2*np.pi, ntheta_total)
    R, THETA = np.meshgrid(r, theta)
    
    # Convert to Cartesian
    X = CYLINDER_CENTER[0] + R * np.cos(THETA)
    Y = CYLINDER_CENTER[1] + R * np.sin(THETA)
    
    # Count point types
    print("\nPoint type distribution:")
    for pt_val, pt_label in POINT_TYPE_LABELS.items():
        count = np.sum(point_types == pt_val)
        percentage = 100.0 * count / point_types.size
        print(f"  {pt_val}: {pt_label:25s} - {count:5d} points ({percentage:5.2f}%)")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Create custom colormap
    cmap = colors.ListedColormap(POINT_TYPE_COLORS)
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    
    # Plot as scatter
    scatter = ax.scatter(X.flatten(), Y.flatten(), c=point_types.flatten(), 
                        cmap=cmap, norm=norm, s=20, marker='o', edgecolors='none')
    
    # Add colorbar
    cbar = fig.colorbar(scatter, ax=ax, ticks=[0, 1, 2, 3, 4, 5])
    cbar.ax.set_yticklabels([POINT_TYPE_LABELS[i] for i in range(6)])
    
    # Draw cylinder (inner boundary - no-slip wall)
    cylinder = plt.Circle(CYLINDER_CENTER, CYLINDER_RADIUS, 
                         fill=False, edgecolor='red', linewidth=3,
                         label=f'Cylinder Surface (no-slip)')
    ax.add_patch(cylinder)
    
    # Draw component outer boundary (interpolation boundary)
    comp_outer = plt.Circle(CYLINDER_CENTER, COMP_OUTER_RADIUS, 
                           fill=False, edgecolor='blue', linewidth=2, linestyle='--',
                           label=f'Outer Boundary (interpolation)')
    ax.add_patch(comp_outer)
    
    # Formatting
    ax.set_xlabel('x (m)', fontsize=12)
    ax.set_ylabel('y (m)', fontsize=12)
    ax.set_title('Circular Component Mesh - Point Types', fontsize=14, weight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=10)
    
    # Set limits
    lim = COMP_OUTER_RADIUS * 1.5
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    
    # Save figure
    out_path = output_dir / "compHoleCut.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\n✓ Saved: {out_path}")
    
    return ntheta_total, nr_total, point_types


def plot_combined_overview(output_dir, comp_dims=None):
    """Plot combined overview showing both meshes together."""
    print("\n=== Combined Overset Mesh Overview ===")
    
    # Load both point type files
    bg_pt = np.loadtxt(output_dir / "bgPtType.txt")
    comp_pt = np.loadtxt(output_dir / "compPtType.txt")
    
    # Background mesh grid
    ny_bg, nx_bg = BG_NY + 2, BG_NX + 2
    bg_pt = bg_pt.reshape(ny_bg, nx_bg)
    dx = (BG_XMAX - BG_XMIN) / BG_NX
    dy = (BG_YMAX - BG_YMIN) / BG_NY
    x_bg = np.linspace(BG_XMIN - dx/2, BG_XMAX + dx/2, nx_bg)
    y_bg = np.linspace(BG_YMIN - dy/2, BG_YMAX + dy/2, ny_bg)
    X_bg, Y_bg = np.meshgrid(x_bg, y_bg)
    
    # Circular mesh grid - use dimensions from previous plot if available
    if comp_dims:
        ntheta_total, nr_total, _ = comp_dims
    else:
        # Infer from data size
        total_points = comp_pt.size
        for nr in range(COMP_NR - 2, COMP_NR + 3):
            for nt in range(COMP_NTHETA - 2, COMP_NTHETA + 3):
                if nr * nt == total_points:
                    ntheta_total, nr_total = nt, nr
                    break
    
    comp_pt = comp_pt.reshape(ntheta_total, nr_total)
    r = np.linspace(CYLINDER_RADIUS, COMP_OUTER_RADIUS, nr_total)
    theta = np.linspace(0, 2*np.pi, ntheta_total)
    R, THETA = np.meshgrid(r, theta)
    X_comp = CYLINDER_CENTER[0] + R * np.cos(THETA)
    Y_comp = CYLINDER_CENTER[1] + R * np.sin(THETA)
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    
    cmap = colors.ListedColormap(POINT_TYPE_COLORS)
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    
    # Left: Background mesh with cylinder region highlighted
    ax1.scatter(X_bg.flatten(), Y_bg.flatten(), c=bg_pt.flatten(), 
               cmap=cmap, norm=norm, s=15, marker='s', edgecolors='none', alpha=0.8)
    cylinder1 = plt.Circle(CYLINDER_CENTER, CYLINDER_RADIUS, 
                          fill=False, edgecolor='black', linewidth=3, linestyle='--')
    ax1.add_patch(cylinder1)
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('y (m)')
    ax1.set_title('Background Rectangular Mesh\n(Hole-cut around cylinder)', weight='bold')
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(-3, 5)
    ax1.set_ylim(-3, 3)
    
    # Right: Circular component mesh
    ax2.scatter(X_comp.flatten(), Y_comp.flatten(), c=comp_pt.flatten(), 
               cmap=cmap, norm=norm, s=30, marker='o', edgecolors='none', alpha=0.8)
    cylinder2 = plt.Circle(CYLINDER_CENTER, CYLINDER_RADIUS, 
                          fill=False, edgecolor='red', linewidth=3)
    comp_outer = plt.Circle(CYLINDER_CENTER, COMP_OUTER_RADIUS, 
                           fill=False, edgecolor='blue', linewidth=2, linestyle='--')
    ax2.add_patch(cylinder2)
    ax2.add_patch(comp_outer)
    ax2.set_xlabel('x (m)')
    ax2.set_ylabel('y (m)')
    ax2.set_title('Circular Component Mesh\n(Around cylinder)', weight='bold')
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    lim = COMP_OUTER_RADIUS * 1.3
    ax2.set_xlim(-lim, lim)
    ax2.set_ylim(-lim, lim)
    
    # Add shared colorbar
    cbar = fig.colorbar(ax1.collections[0], ax=[ax1, ax2], 
                       ticks=[0, 1, 2, 3, 4, 5], fraction=0.02)
    cbar.ax.set_yticklabels([POINT_TYPE_LABELS[i] for i in range(6)])
    
    fig.suptitle('Overset Mesh Configuration - Cylinder Flow (Re=100)', 
                fontsize=16, weight='bold', y=0.98)
    
    # Save
    out_path = output_dir / "overset_mesh_overview.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\n✓ Saved: {out_path}")


def main():
    # Find cylinder flow directory
    output_dirs = [d for d in TEMP_DIR.iterdir() 
                  if d.is_dir() and d.name.startswith("CylinderFlow")]
    
    if not output_dirs:
        print("Error: No CylinderFlow output directory found in TEMP/")
        print("Please run overset_cylinder.exe first!")
        return
    
    # Use the most recently modified directory
    output_dir = max(output_dirs, key=lambda d: d.stat().st_mtime)
    print(f"Reading hole cutting data from: {output_dir.name}\n")
    
    # Generate all visualizations
    plot_background_hole_cutting(output_dir)
    comp_dims = plot_circular_mesh_connectivity(output_dir)
    if comp_dims and comp_dims[0] is not None:
        plot_combined_overview(output_dir, comp_dims)
    else:
        print("\nSkipping combined overview due to circular mesh dimension issues")
    
    print("\n" + "="*60)
    print("✓ All hole cutting visualizations generated successfully!")
    print(f"  Location: {output_dir}")
    print("="*60)


if __name__ == "__main__":
    main()
