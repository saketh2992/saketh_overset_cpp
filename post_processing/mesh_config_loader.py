"""
Mesh configuration loader for Python post-processing scripts.
Reads mesh_config.json to ensure consistency with C++ solver.
"""
from __future__ import annotations
import json
import numpy as np
from pathlib import Path
from typing import Tuple, List

DTYPE = np.float64
PI = np.pi


class MeshConfig:
    """Load and provide access to mesh configuration from JSON file."""
    
    def __init__(self, config_file: Path | str = None):
        """
        Load mesh configuration from JSON file.
        
        Args:
            config_file: Path to mesh_config.json. If None, looks in workspace root.
        """
        if config_file is None:
            # Default to workspace root
            script_dir = Path(__file__).resolve().parent
            config_file = script_dir.parent / "mesh_config.json"
        else:
            config_file = Path(config_file)
            
        if not config_file.exists():
            raise FileNotFoundError(f"Mesh configuration file not found: {config_file}")
        
        with open(config_file, 'r') as f:
            self._config = json.load(f)
        
        self.mulFac = self._config['mulFac']
        self.meshes = self._config['meshes']
        self.solver = self._config.get('solver', {})
    
    def get_mesh_desc(self, mesh_name: str) -> Tuple[float, float, float, float, float, int, int]:
        """
        Get mesh descriptor tuple for a specific mesh.
        
        Args:
            mesh_name: Name of mesh ('bg' or 'comp')
            
        Returns:
            Tuple: (x0, y0, theta_rad, length, width, Nx, Ny)
        """
        mesh = self.meshes[mesh_name]
        return (
            mesh['x0'],
            mesh['y0'],
            mesh['theta_deg'] * PI / 180.0,  # Convert to radians
            mesh['length'],
            mesh['width'],
            mesh['Nx_base'] * self.mulFac,
            mesh['Ny_base'] * self.mulFac
        )
    
    def get_all_mesh_descriptors(self) -> List[Tuple[float, float, float, float, float, int, int]]:
        """
        Get mesh descriptors for all meshes in order: [bg, comp].
        
        Returns:
            List of tuples: [(x0, y0, theta_rad, length, width, Nx, Ny), ...]
        """
        return [
            self.get_mesh_desc('bg'),
            self.get_mesh_desc('comp')
        ]
    
    def get_mesh_names(self) -> List[str]:
        """Get list of mesh names in order."""
        return ['bg', 'comp']
    
    def get_var_names(self) -> List[str]:
        """Get list of variable names."""
        return self.solver.get('variables', ['U', 'V', 'P'])
    
    def __repr__(self):
        return f"MeshConfig(mulFac={self.mulFac}, meshes={list(self.meshes.keys())})"


# Convenience function for quick access
def load_mesh_config(config_file: Path | str = None) -> MeshConfig:
    """
    Load mesh configuration.
    
    Args:
        config_file: Path to mesh_config.json. If None, looks in workspace root.
        
    Returns:
        MeshConfig object
    """
    return MeshConfig(config_file)


if __name__ == "__main__":
    # Test loading
    config = load_mesh_config()
    print(config)
    print(f"mulFac: {config.mulFac}")
    print(f"BG mesh: {config.get_mesh_desc('bg')}")
    print(f"Comp mesh: {config.get_mesh_desc('comp')}")
    print(f"All descriptors: {config.get_all_mesh_descriptors()}")
