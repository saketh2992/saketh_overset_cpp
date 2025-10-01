#!/bin/bash
# Script to download, build, and install TIOGA library

set -e  # Exit on error

echo "==================================="
echo "TIOGA Library Build Script"
echo "==================================="

# Configuration
INSTALL_DIR="${HOME}/tioga_install"
BUILD_DIR="./tioga_build_temp"

echo ""
echo "Installation directory: ${INSTALL_DIR}"
echo "Build directory: ${BUILD_DIR}"
echo ""

# Ensure cmake is in PATH
export PATH="/usr/local/bin:/opt/homebrew/bin:$PATH"

# Check for required tools
echo "Checking for required tools..."
command -v git >/dev/null 2>&1 || { echo "Error: git is required but not installed."; exit 1; }
command -v cmake >/dev/null 2>&1 || { echo "Error: cmake is required but not installed."; exit 1; }
echo "✓ All required tools found"
echo ""

# Check for MPI (required by TIOGA)
echo "Checking for MPI..."
if command -v mpicc >/dev/null 2>&1; then
    echo "✓ MPI found: $(mpicc --version | head -n 1)"
else
    echo "⚠ Warning: MPI not found. Installing via Homebrew..."
    brew install open-mpi
fi
echo ""

# Clean up old build directory if it exists
if [ -d "${BUILD_DIR}" ]; then
    echo "Removing old build directory..."
    rm -rf "${BUILD_DIR}"
fi

# Clone TIOGA repository
echo "Cloning TIOGA repository..."
git clone https://github.com/jsitaraman/tioga.git "${BUILD_DIR}"
cd "${BUILD_DIR}"
echo "✓ TIOGA cloned successfully"
echo ""

# Create build directory
mkdir -p build
cd build

# Configure with CMake
echo "Configuring TIOGA with CMake..."
cmake .. \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=11 \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpicxx

echo "✓ Configuration complete"
echo ""

# Build TIOGA
echo "Building TIOGA (this may take a few minutes)..."
make -j$(sysctl -n hw.ncpu)
echo "✓ Build complete"
echo ""

# Install TIOGA
echo "Installing TIOGA to ${INSTALL_DIR}..."
make install
echo "✓ Installation complete"
echo ""

# Go back to original directory
cd ../..

# Clean up build directory
echo "Cleaning up temporary build directory..."
rm -rf "${BUILD_DIR}"
echo "✓ Cleanup complete"
echo ""

echo "==================================="
echo "TIOGA Installation Summary"
echo "==================================="
echo "Installation directory: ${INSTALL_DIR}"
echo ""
echo "Library files:"
ls -lh "${INSTALL_DIR}/lib/"
echo ""
echo "Include files:"
ls -lh "${INSTALL_DIR}/include/"
echo ""
echo "==================================="
echo "✓ TIOGA installation successful!"
echo "==================================="
echo ""
echo "Next steps:"
echo "1. Update CMakeLists.txt with TIOGA configuration"
echo "2. Build your solver with: cmake --preset=release && cmake --build build/release"
echo "3. TIOGA will be linked automatically"
echo ""
