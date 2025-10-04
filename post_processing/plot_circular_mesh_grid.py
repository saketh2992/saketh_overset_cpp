"""
Detailed visualization of circular component mesh showing grid structure.
- Shows all grid points clearly with lines connecting them
- Displays both radial and circumferential grid lines
- Color codes points by type
- Useful for understanding the polar mesh structure

Usage (PowerShell):
  py post_processing/plot_circular_mesh_grid.py
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors
import numpy as np
from pathlib import Path

# Configuration
ROOT = Path(__file__).resolve().parents[1]
TEMP_DIR = ROOT / "TEMP"

# Circular mesh parameters (from main_cylinder_flow.cpp)
CYLINDER_CENTER = (0.0, 0.0)
CYLINDER_RADIUS = 0.5
COMP_OUTER_RADIUS = 2.0  # Updated outer radius
COMP_NTHETA = 25
COMP_NR = 25

# Point type labels
POINT_TYPE_LABELS = {
    0: "UNUSED",
    1: "CALCULATED",
    2: "DONOR",
    3: "INTERPOLATION_RECIEVER",
    4: "DONOR_BUFFER",
    5: "BC_SPECIFIED"
}

POINT_TYPE_COLORS = ['gray', 'lightblue', 'orange', 'red', 'yellow', 'green']

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 10})


def plot_circular_mesh_detailed(output_dir):
    """Create detailed visualization of circular mesh with grid lines."""
    print("\n=== Detailed Circular Mesh Grid Visualization ===")
    
    # Load point type data
    pt_file = output_dir / "compPtType.txt"
    if not pt_file.exists():
        print(f"Error: {pt_file} not found!")
        return
    
    point_types = np.loadtxt(pt_file)
    print(f"Point type data size: {point_types.size}")
    
    # Infer grid dimensions
    total_points = point_types.size
    possible_dims = []
    for nr in range(COMP_NR - 2, COMP_NR + 3):
        for nt in range(COMP_NTHETA - 2, COMP_NTHETA + 3):
            if nr * nt == total_points:
                possible_dims.append((nt, nr))
    
    if not possible_dims:
        print(f"Error: Cannot determine grid dimensions for {total_points} points")
        return
    
    ntheta_total, nr_total = possible_dims[0]
    print(f"Using grid dimensions: {ntheta_total} (theta) x {nr_total} (radial)")
    point_types = point_types.reshape(ntheta_total, nr_total)
    
    # Generate circular mesh coordinates
    r = np.linspace(CYLINDER_RADIUS, COMP_OUTER_RADIUS, nr_total)
    theta = np.linspace(0, 2*np.pi, ntheta_total, endpoint=False)
    
    # Create meshgrid
    R, THETA = np.meshgrid(r, theta)
    X = CYLINDER_CENTER[0] + R * np.cos(THETA)
    Y = CYLINDER_CENTER[1] + R * np.sin(THETA)
    
    # Create figure with larger size for clarity
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9))
    
    # --- LEFT PLOT: Grid structure with lines ---
    print("\nPlotting grid structure...")
    
    # Plot radial lines (constant theta)
    for i in range(0, ntheta_total, max(1, ntheta_total // 12)):  # Show ~12 radial lines
        ax1.plot(X[i, :], Y[i, :], 'b-', linewidth=0.5, alpha=0.3)
    
    # Plot circumferential lines (constant r)
    for j in range(nr_total):
        # Close the circle by adding first point at the end
        x_circle = np.append(X[:, j], X[0, j])
        y_circle = np.append(Y[:, j], Y[0, j])
        ax1.plot(x_circle, y_circle, 'b-', linewidth=0.5, alpha=0.3)
    
    # Plot all grid points
    cmap = colors.ListedColormap(POINT_TYPE_COLORS)
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    
    scatter1 = ax1.scatter(X.flatten(), Y.flatten(), c=point_types.flatten(), 
                          cmap=cmap, norm=norm, s=15, marker='o', 
                          edgecolors='black', linewidths=0.3, zorder=5)
    
    # Draw cylinder (inner boundary)
    cylinder1 = plt.Circle(CYLINDER_CENTER, CYLINDER_RADIUS, 
                          fill=False, edgecolor='red', linewidth=3,
                          linestyle='--', label='Cylinder surface (no-slip)')
    ax1.add_patch(cylinder1)
    
    # Draw outer boundary
    outer1 = plt.Circle(CYLINDER_CENTER, COMP_OUTER_RADIUS, 
                       fill=False, edgecolor='blue', linewidth=2,
                       linestyle='--', label='Outer boundary (interpolation)')
    ax1.add_patch(outer1)
    
    # Add grid info
    ax1.text(0.02, 0.98, f'Grid: {ntheta_total}θ × {nr_total}r = {total_points} points\n'
                         f'Δr = {(COMP_OUTER_RADIUS - CYLINDER_RADIUS) / (nr_total - 1):.4f}\n'
                         f'Δθ = {360.0 / ntheta_total:.2f}°',
             transform=ax1.transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax1.set_xlabel('x (m)', fontsize=12)
    ax1.set_ylabel('y (m)', fontsize=12)
    ax1.set_title('Circular Mesh - Grid Structure with Lines', fontsize=14, weight='bold')
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.2, linestyle=':')
    ax1.legend(loc='upper right', fontsize=9)
    
    # Set limits
    lim = COMP_OUTER_RADIUS * 1.15
    ax1.set_xlim(-lim, lim)
    ax1.set_ylim(-lim, lim)
    
    # --- RIGHT PLOT: Zoomed view of a section ---
    print("Plotting zoomed section...")
    
    # Show a 90-degree sector zoomed in
    sector_indices = slice(0, max(1, ntheta_total // 4))
    
    # Plot grid lines in sector
    for i in range(0, ntheta_total // 4, max(1, ntheta_total // 24)):
        ax2.plot(X[i, :], Y[i, :], 'b-', linewidth=1, alpha=0.5)
    
    for j in range(nr_total):
        x_arc = X[sector_indices, j]
        y_arc = Y[sector_indices, j]
        ax2.plot(x_arc, y_arc, 'b-', linewidth=1, alpha=0.5)
    
    # Plot points in sector with labels
    X_sector = X[sector_indices, :]
    Y_sector = Y[sector_indices, :]
    pt_sector = point_types[sector_indices, :]
    
    scatter2 = ax2.scatter(X_sector.flatten(), Y_sector.flatten(), 
                          c=pt_sector.flatten(), cmap=cmap, norm=norm, 
                          s=50, marker='o', edgecolors='black', linewidths=0.5, zorder=5)
    
    # Add point numbers for some points
    step = max(1, nr_total // 5)
    for i in range(0, X_sector.shape[0], max(1, X_sector.shape[0] // 6)):
        for j in range(0, nr_total, step):
            ax2.text(X_sector[i, j], Y_sector[i, j], f'{i},{j}', 
                    fontsize=6, ha='center', va='bottom', alpha=0.7)
    
    # Draw boundaries
    cylinder2 = plt.Circle(CYLINDER_CENTER, CYLINDER_RADIUS, 
                          fill=False, edgecolor='red', linewidth=2, linestyle='--')
    ax2.add_patch(cylinder2)
    
    outer2 = plt.Circle(CYLINDER_CENTER, COMP_OUTER_RADIUS, 
                       fill=False, edgecolor='blue', linewidth=2, linestyle='--')
    ax2.add_patch(outer2)
    
    # Draw sector lines
    theta_max = 2 * np.pi / 4
    r_line = np.linspace(0, COMP_OUTER_RADIUS * 1.1, 10)
    ax2.plot([0, r_line[-1]], [0, 0], 'k--', linewidth=1, alpha=0.5)
    ax2.plot(r_line * np.cos(theta_max), r_line * np.sin(theta_max), 
            'k--', linewidth=1, alpha=0.5)
    
    ax2.set_xlabel('x (m)', fontsize=12)
    ax2.set_ylabel('y (m)', fontsize=12)
    ax2.set_title('Zoomed View - First Quadrant\n(showing grid indices)', 
                 fontsize=14, weight='bold')
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    
    ax2.set_xlim(-0.2, COMP_OUTER_RADIUS * 1.1)
    ax2.set_ylim(-0.2, COMP_OUTER_RADIUS * 1.1)
    
    # Add shared colorbar
    cbar = fig.colorbar(scatter1, ax=[ax1, ax2], ticks=[0, 1, 2, 3, 4, 5], 
                       fraction=0.02, pad=0.04)
    cbar.ax.set_yticklabels([POINT_TYPE_LABELS[i] for i in range(6)])
    cbar.set_label('Point Type', fontsize=11)
    
    # Overall title
    fig.suptitle(f'Circular Component Mesh Detail - Re=100\n'
                f'Inner R={CYLINDER_RADIUS}, Outer R={COMP_OUTER_RADIUS}, '
                f'{ntheta_total}θ × {nr_total}r points', 
                fontsize=16, weight='bold', y=0.98)
    
    # Count point types
    print("\nPoint type distribution:")
    for pt_val, pt_label in POINT_TYPE_LABELS.items():
        count = np.sum(point_types == pt_val)
        percentage = 100.0 * count / point_types.size
        print(f"  {pt_val}: {pt_label:25s} - {count:5d} points ({percentage:5.2f}%)")
    
    # Save
    out_path = output_dir / "circular_mesh_detail.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=250, bbox_inches='tight')
    plt.close(fig)
    print(f"\n✓ Saved: {out_path}")


def plot_mesh_3d_view(output_dir):
    """Create a 3D-style view showing the mesh structure."""
    print("\n=== Creating 3D perspective view ===")
    
    # Load point type data
    pt_file = output_dir / "compPtType.txt"
    if not pt_file.exists():
        return
    
    point_types = np.loadtxt(pt_file)
    total_points = point_types.size
    
    # Infer dimensions
    for nr in range(COMP_NR - 2, COMP_NR + 3):
        for nt in range(COMP_NTHETA - 2, COMP_NTHETA + 3):
            if nr * nt == total_points:
                ntheta_total, nr_total = nt, nr
                break
    
    point_types = point_types.reshape(ntheta_total, nr_total)
    
    # Generate mesh
    r = np.linspace(CYLINDER_RADIUS, COMP_OUTER_RADIUS, nr_total)
    theta = np.linspace(0, 2*np.pi, ntheta_total, endpoint=False)
    R, THETA = np.meshgrid(r, theta)
    X = CYLINDER_CENTER[0] + R * np.cos(THETA)
    Y = CYLINDER_CENTER[1] + R * np.sin(THETA)
    
    # Create figure
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111)
    
    # Plot mesh with different colors for boundaries
    cmap = colors.ListedColormap(POINT_TYPE_COLORS)
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    
    # Plot all radial lines
    for i in range(ntheta_total):
        ax.plot(X[i, :], Y[i, :], 'b-', linewidth=0.8, alpha=0.4, zorder=1)
    
    # Plot all circumferential lines
    for j in range(nr_total):
        x_circle = np.append(X[:, j], X[0, j])
        y_circle = np.append(Y[:, j], Y[0, j])
        ax.plot(x_circle, y_circle, 'b-', linewidth=0.8, alpha=0.4, zorder=1)
    
    # Plot points
    scatter = ax.scatter(X.flatten(), Y.flatten(), c=point_types.flatten(), 
                        cmap=cmap, norm=norm, s=25, marker='o', 
                        edgecolors='black', linewidths=0.5, zorder=3)
    
    # Highlight boundaries
    # Inner boundary (cylinder)
    ax.plot(X[:, 0], Y[:, 0], 'ro', markersize=6, label='Inner (no-slip)', zorder=4)
    # Outer boundary (interpolation)
    ax.plot(X[:, -1], Y[:, -1], 'bs', markersize=6, label='Outer (interpolation)', zorder=4)
    
    # Draw cylinder
    cylinder = plt.Circle(CYLINDER_CENTER, CYLINDER_RADIUS, 
                         fill=True, facecolor='white', edgecolor='red', 
                         linewidth=3, zorder=5, alpha=0.9)
    ax.add_patch(cylinder)
    
    # Add text
    ax.text(0, 0, 'Cylinder\n(No-slip wall)', ha='center', va='center', 
           fontsize=11, weight='bold', zorder=6)
    
    # Colorbar
    cbar = fig.colorbar(scatter, ax=ax, ticks=[0, 1, 2, 3, 4, 5], fraction=0.046)
    cbar.ax.set_yticklabels([POINT_TYPE_LABELS[i] for i in range(6)])
    
    ax.set_xlabel('x (m)', fontsize=13)
    ax.set_ylabel('y (m)', fontsize=13)
    ax.set_title(f'Circular Mesh - Complete Grid Structure\n'
                f'{ntheta_total} circumferential × {nr_total} radial = {total_points} points',
                fontsize=15, weight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(loc='upper right', fontsize=11)
    
    lim = COMP_OUTER_RADIUS * 1.1
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    
    # Save
    out_path = output_dir / "circular_mesh_complete.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=250, bbox_inches='tight')
    plt.close(fig)
    print(f"✓ Saved: {out_path}")


def main():
    # Find latest cylinder flow directory
    output_dirs = [d for d in TEMP_DIR.iterdir() 
                  if d.is_dir() and d.name.startswith("CylinderFlow")]
    
    if not output_dirs:
        print("Error: No CylinderFlow output directory found in TEMP/")
        print("Please run overset_cylinder.exe first!")
        return
    
    output_dir = max(output_dirs, key=lambda d: d.stat().st_mtime)
    print(f"Reading mesh data from: {output_dir.name}\n")
    
    # Generate visualizations
    plot_circular_mesh_detailed(output_dir)
    plot_mesh_3d_view(output_dir)
    
    print("\n" + "="*70)
    print("✓ Circular mesh grid visualizations complete!")
    print(f"  Location: {output_dir}")
    print("  Files: circular_mesh_detail.png, circular_mesh_complete.png")
    print("="*70)


if __name__ == "__main__":
    main()
