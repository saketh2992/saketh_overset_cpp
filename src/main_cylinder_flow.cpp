#include <iostream>
#include <sstream>
#include <iomanip>
#ifdef _WIN32
#include <direct.h>  // For _mkdir on Windows
#else
#include <sys/stat.h>  // For mkdir on Unix/Linux/macOS
#endif
#include "util/constants.h"
#include "util/utilities.h"
#include "mesh/adt.h"
#include "mesh/datastructure.h"
#include "mesh/datastructure_rectangular.h"
#include "mesh/datastructure_circular.h"
#include "solver/solver.h"
#include "io/output.h"

using namespace std;

// Helper function to create directory (cross-platform)
void createDirectory(const string& path) {
#ifdef _WIN32
    _mkdir(path.c_str());
#else
    mkdir(path.c_str(), 0755);
#endif
}

int main() {
    cout << "===========================================================" << endl;
    cout << "    FLOW PAST CIRCULAR CYLINDER - OVERSET MESH SOLVER     " << endl;
    cout << "    Background: Rectangular | Component: Circular         " << endl;
    cout << "===========================================================" << endl;
    
    // ============================================================
    // SIMULATION PARAMETERS
    // ============================================================
    const double Re = 40.0;              // Reynolds number
    const double U_inf = 1.0;             // Freestream velocity
    
    // Background rectangular mesh parameters
    const double bg_x0 = -5.0;            // Bottom-left corner X
    const double bg_y0 = -5.0;            // Bottom-left corner Y
    const double bg_Lx = 15.0;            // Length in X direction (-5 to 10)
    const double bg_Ly = 10.0;            // Length in Y direction (-5 to 5)
    const int bg_Nx = 50;                 // Grid points in X
    const int bg_Ny = 50;                 // Grid points in Y
    const double bg_theta = 0.0;          // No rotation
    
    // Circular mesh parameters (around cylinder)
    const double cyl_centerX = 0.0;       // Cylinder center X
    const double cyl_centerY = 0.0;       // Cylinder center Y
    const double cyl_radius = 0.5;        // Cylinder radius
    const double comp_outerRadius = 2.0;  // Component mesh outer radius
    const int comp_Ntheta = 25;           // Points in circumferential direction
    const int comp_Nr = 25;               // Points in radial direction
    
    const int maxIterations = 50000;      // Maximum solver iterations
    
    cout << "\n--------------------------------------------------" << endl;
    cout << "Simulation Parameters:" << endl;
    cout << "--------------------------------------------------" << endl;
    cout << "Reynolds Number: " << Re << endl;
    cout << "Freestream Velocity: " << U_inf << " m/s" << endl;
    cout << "\nBackground Mesh (Rectangular):" << endl;
    cout << "  Domain: [" << bg_x0 << ", " << bg_x0 + bg_Lx << "] x [" 
         << bg_y0 << ", " << bg_y0 + bg_Ly << "]" << endl;
    cout << "  Grid: " << bg_Nx << " x " << bg_Ny << " points" << endl;
    cout << "  dx = " << bg_Lx / bg_Nx << ", dy = " << bg_Ly / bg_Ny << endl;
    cout << "\nCylinder & Component Mesh (Circular):" << endl;
    cout << "  Center: (" << cyl_centerX << ", " << cyl_centerY << ")" << endl;
    cout << "  Cylinder Radius: " << cyl_radius << endl;
    cout << "  Component Outer Radius: " << comp_outerRadius << endl;
    cout << "  Grid: " << comp_Ntheta << " (theta) x " << comp_Nr << " (radial) points" << endl;
    cout << "  dtheta = " << 360.0 / comp_Ntheta << " degrees" << endl;
    cout << "  dr = " << (comp_outerRadius - cyl_radius) / comp_Nr << endl;
    cout << "--------------------------------------------------" << endl;
    
    // ============================================================
    // CREATE MESHES
    // ============================================================
    cout << "\n--------------------------------------------------" << endl;
    cout << "Creating meshes..." << endl;
    cout << "--------------------------------------------------" << endl;
    
    // Background rectangular mesh
    cout << "Creating background rectangular mesh..." << endl;
    DataStructureRectangular bgMesh(bg_x0, bg_y0, bg_theta,
                                     bg_Lx, bg_Ly, 
                                     bg_Nx, bg_Ny, Re);
    bgMesh.flowType = "CYLINDER_FLOW";  // Set flow type
    bgMesh.uInlet = U_inf;              // Set inlet velocity
    
    // Component circular mesh around cylinder
    cout << "Creating circular mesh around cylinder..." << endl;
    DataStructureCircular compMesh(cyl_centerX, cyl_centerY,
                                    cyl_radius, comp_outerRadius,
                                    comp_Nr, comp_Ntheta, Re);
    compMesh.flowType = "CYLINDER_FLOW";
    
    // Create output directory
    ostringstream dirName;
    dirName << "TEMP/CylinderFlow_Re" << static_cast<int>(Re)
            << "_bg" << bg_Nx << "x" << bg_Ny
            << "_cyl" << comp_Ntheta << "x" << comp_Nr
            << "_R" << fixed << setprecision(2) << cyl_radius;
    
    string outputDir = dirName.str();
    createDirectory("TEMP");
    createDirectory(outputDir);
    
    cout << "Output directory: " << outputDir << endl;
    
    // ============================================================
    // OVERSET MESH CONNECTIVITY
    // ============================================================
    cout << "\n--------------------------------------------------" << endl;
    cout << "Setting up overset mesh connectivity..." << endl;
    cout << "--------------------------------------------------" << endl;
    
    cout << "1. Marking boundary donors on circular mesh..." << endl;
    compMesh.MarkBoundaryDonor(bgMesh);
    
    cout << "2. Performing hole cutting (removing points inside cylinder)..." << endl;
    bgMesh.HoleCutting(compMesh);
    
    cout << "3. Marking boundary point types..." << endl;
    // Background mesh: Inlet, outlet, and far-field boundaries
    bgMesh.MarkBoundaryPointType(BC_SPECIFIED);
    // Circular mesh outer boundary: receives interpolation from background
    // Circular mesh inner boundary (cylinder surface): no-slip wall
    compMesh.MarkBoundaryPointType(INTERPOLATION_RECIEVER);
    
    cout << "4. Storing internal interpolated points..." << endl;
    bgMesh.StoreInternalInterpolatedPoints();
    compMesh.StoreInternalInterpolatedPoints();
    
    cout << "5. Writing point type files..." << endl;
    bgMesh.WriteTxtPointType(outputDir + "/bgPtType.txt");
    compMesh.WriteTxtPointType(outputDir + "/compPtType.txt");
    
    // ============================================================
    // INITIALIZE AND SOLVE
    // ============================================================
    cout << "\n--------------------------------------------------" << endl;
    cout << "Initializing flow field..." << endl;
    cout << "--------------------------------------------------" << endl;
    
    // Initialize with freestream conditions
    cout << "Initializing background mesh (freestream: U=" << U_inf << ", V=0, P=0)..." << endl;
    initialize(&bgMesh, &compMesh);
    
    cout << "Initializing circular mesh (freestream conditions)..." << endl;
    initialize(&compMesh, &bgMesh);
    
    cout << "\n--------------------------------------------------" << endl;
    cout << "Solving Navier-Stokes equations..." << endl;
    cout << "--------------------------------------------------" << endl;
    
    // Solve using cylinder-specific boundary conditions
    SolveCylinder(&bgMesh, &compMesh, maxIterations, U_inf);
    
    // ============================================================
    // POST-PROCESSING
    // ============================================================
    cout << "\n--------------------------------------------------" << endl;
    cout << "Writing output files..." << endl;
    cout << "--------------------------------------------------" << endl;
    
    GetOutput(&bgMesh, (outputDir + "/output_bgMesh.dat").c_str());
    GetOutput(&compMesh, (outputDir + "/output_circMesh.dat").c_str());
    
    // Write visualization files
    // NOTE: WriteGnuplotAll expects rectangular meshes, skip for now with circular mesh
    // WriteGnuplotAll(bgMesh, compMesh, 
    //                 "Flow Past Circular Cylinder Re=" + to_string((int)Re), 
    //                 "cylinder_flow", outputDir);
    cout << "Basic output files written (detailed plots skipped for circular mesh)" << endl;
    
    // ============================================================
    // CALCULATE AERODYNAMIC FORCES
    // ============================================================
    cout << "\n--------------------------------------------------" << endl;
    cout << "Calculating aerodynamic forces..." << endl;
    cout << "--------------------------------------------------" << endl;
    
    // Calculate drag and lift coefficients
    // This would require surface integration on the cylinder
    // For now, just placeholder
    cout << "Drag Coefficient (Cd): [To be implemented]" << endl;
    cout << "Lift Coefficient (Cl): [To be implemented]" << endl;
    cout << "Strouhal Number (St): [To be implemented]" << endl;
    
    cout << "\n===========================================================" << endl;
    cout << "    SIMULATION COMPLETE!                                    " << endl;
    cout << "    Results saved in: " << outputDir << endl;
    cout << "===========================================================" << endl;
    cout << "\nExpected flow features at Re=" << Re << ":" << endl;
    cout << "  - Symmetric steady flow (if Re < 40)" << endl;
    cout << "  - Vortex pair in wake (if 40 < Re < 47)" << endl;
    cout << "  - Vortex shedding (if Re > 47, Karman vortex street)" << endl;
    cout << "===========================================================" << endl;
    
    return 0;
}
