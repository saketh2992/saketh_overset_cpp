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
#include "util/mesh_config.h"
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
    cout << "==================================================" << endl;
    cout << "    OVERSET MESH SOLVER - NEW IMPLEMENTATION     " << endl;
    cout << "    Using Inheritance-Based DataStructure        " << endl;
    cout << "==================================================" << endl;
    
    // Load mesh configuration from JSON file
    MeshConfig config("mesh_config.json");
    
    // Create meshes using the new inheritance-based structure
    // Using DataStructureRectangular for both background and component meshes
    DataStructureRectangular bgMesh(config.getBgX0(), config.getBgY0(), config.getBgTheta(),
                                     config.getBgLength(), config.getBgWidth(), 
                                     config.getBgNx(), config.getBgNy(), config.getReynolds());
    
    DataStructureRectangular compMesh(config.getCompX0(), config.getCompY0(), config.getCompTheta(),
                                       config.getCompLength(), config.getCompWidth(),
                                       config.getCompNx(), config.getCompNy(), config.getReynolds());

    // Create output directory name based on parameters
    // Format: TEMP/Re{Re}_bg{Nx}x{Ny}_comp{Nx}x{Ny}_theta{angle}
    ostringstream dirName;
    dirName << "TEMP/Re" << static_cast<int>(config.getReynolds())
            << "_bg" << config.getBgNx() << "x" << config.getBgNy()
            << "_comp" << config.getCompNx() << "x" << config.getCompNy()
            << "_theta" << fixed << setprecision(0) << (config.getCompTheta() * 180.0 / PI);
    
    string outputDir = dirName.str();
    
    // Create directories
    createDirectory("TEMP");
    createDirectory(outputDir);
    
    cout << "\nOutput directory: " << outputDir << endl;
    cout << "\n--------------------------------------------------" << endl;
    cout << "Starting overset mesh processing..." << endl;
    cout << "--------------------------------------------------" << endl;

    // Overset mesh connectivity operations
    cout << "\n1. Marking boundary donors..." << endl;
    compMesh.MarkBoundaryDonor(bgMesh);
    
    cout << "2. Performing hole cutting..." << endl;
    bgMesh.HoleCutting(compMesh);
    
    cout << "3. Marking boundary point types..." << endl;
    bgMesh.MarkBoundaryPointType(BC_SPECIFIED);
    compMesh.MarkBoundaryPointType(INTERPOLATION_RECIEVER);
    
    cout << "4. Storing internal interpolated points..." << endl;
    bgMesh.StoreInternalInterpolatedPoints();
    compMesh.StoreInternalInterpolatedPoints();

    cout << "5. Writing point type files..." << endl;
    bgMesh.WriteTxtPointType(outputDir + "/bgPtType.txt");
    compMesh.WriteTxtPointType(outputDir + "/compPtType.txt");

    cout << "\n--------------------------------------------------" << endl;
    cout << "Initializing and solving..." << endl;
    cout << "--------------------------------------------------" << endl;
    
    // Initialize and solve
    cout << "Initializing background mesh..." << endl;
    initialize(&bgMesh, &compMesh);
    
    cout << "Initializing component mesh..." << endl;
    initialize(&compMesh, &bgMesh);
    
    cout << "\nSolving coupled overset mesh system..." << endl;
    Solve(&bgMesh, &compMesh, config.getMaxIterations());
    
    cout << "\n--------------------------------------------------" << endl;
    cout << "Writing output files..." << endl;
    cout << "--------------------------------------------------" << endl;
    
    // Write output
    GetOutput(&bgMesh, (outputDir + "/output_Upwind_bgMesh.dat").c_str());
    GetOutput(&compMesh, (outputDir + "/output_Upwind_compMesh.dat").c_str());
    
    // Also emit gnuplot data and a plotting script
    WriteGnuplotAll(bgMesh, compMesh, "NS contour for 2d plate", "lid_overset", outputDir);
    
    cout << "\n==================================================" << endl;
    cout << "    SIMULATION COMPLETE!                          " << endl;
    cout << "    Results saved in: " << outputDir << endl;
    cout << "==================================================" << endl;

    return 0;
}
