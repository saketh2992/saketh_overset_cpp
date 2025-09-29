# Overset Mesh Solver (C++)

A simple overset mesh CFD solver written in C++.

## Build (MinGW + CMake)

```powershell
# Ensure MinGW is on PATH
$env:Path = "C:\\MinGW\\bin;" + $env:Path

# Configure & build Release
cmake --preset mingw-release
cmake --build --preset mingw-release -j

# Run
.\build\release\overset.exe
```

## Build (one-shot g++)

```powershell
$env:Path = "C:\\MinGW\\bin;" + $env:Path

g++ -std=c++11 -O3 main.cpp adt.cpp datastructure.cpp solver.cpp utilities.cpp output.cpp -o build\overset.exe
.\build\overset.exe
```

## Notes
- Requires MinGW (g++, gcc, mingw32-make).
- CMake presets are provided in `CMakePresets.json`.
- Source code is C++11-compatible (works with GCC 6.x).
