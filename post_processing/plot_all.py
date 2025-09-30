"""
Master post-processing script - runs all visualization scripts.
Generates all plots from solver output in one go:
  - Contour plots (U, V, P)
  - Hole cutting visualization
  - Centerline comparisons with Ghia et al.

Usage (PowerShell):
  py post_processing/plot_all.py [optional: output_directory]
  
Examples:
  py post_processing/plot_all.py                              # Use latest run
  py post_processing/plot_all.py TEMP/Re100_bg50x50_comp30x30_theta0
"""
from __future__ import annotations
import sys
from pathlib import Path

# Add post_processing directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent))

def find_latest_output_dir():
    """Find the most recently modified output directory in TEMP."""
    root = Path(__file__).resolve().parents[1]
    temp_dir = root / "TEMP"
    
    if not temp_dir.exists():
        raise FileNotFoundError("TEMP directory not found. Run the solver first!")
    
    # Find all subdirectories that match the pattern Re*
    output_dirs = [d for d in temp_dir.iterdir() if d.is_dir() and d.name.startswith("Re")]
    
    if not output_dirs:
        raise FileNotFoundError("No output directories found in TEMP/. Run the solver first!")
    
    # Return the most recently modified directory
    latest_dir = max(output_dirs, key=lambda d: d.stat().st_mtime)
    return latest_dir

# Determine output directory
if len(sys.argv) > 1:
    output_dir = Path(sys.argv[1])
    if not output_dir.exists():
        print(f"Error: Directory '{output_dir}' does not exist.")
        sys.exit(1)
else:
    output_dir = find_latest_output_dir()

print("=" * 70)
print("OVERSET MESH POST-PROCESSING - GENERATING ALL PLOTS")
print("=" * 70)
print(f"Using output directory: {output_dir}")
print()

# Import and run each plotting module
try:
    print("1/3 - Generating contour plots...")
    print("-" * 70)
    import plot_contours
    plot_contours.main(output_dir)
    print()
    
    print("2/3 - Generating hole cutting visualizations...")
    print("-" * 70)
    import plot_hole_cutting
    plot_hole_cutting.main(output_dir)
    print()
    
    print("3/3 - Generating centerline plots...")
    print("-" * 70)
    import plot_centerline
    plot_centerline.main(output_dir)
    print()
    
    print("=" * 70)
    print("✓ ALL PLOTS GENERATED SUCCESSFULLY!")
    print("=" * 70)
    print()
    print(f"Output location: {output_dir}/")
    print("Generated files:")
    print("  • Contours: lid_overset_U.png, lid_overset_V.png, lid_overset_P.png")
    print("  • Hole cutting: bgHoleCut.png, compHoleCut.png")
    print("  • Centerlines: lid_overset_centerline_U.png, lid_overset_centerline_V.png")
    
except Exception as e:
    print()
    print("=" * 70)
    print("✗ ERROR DURING POST-PROCESSING")
    print("=" * 70)
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
