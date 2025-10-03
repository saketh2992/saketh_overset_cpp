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
    // Load mesh configuration from JSON file
    MeshConfig config("mesh_config.json");
    
    DataStructure bgMesh(config.getBgX0(), config.getBgY0(), config.getBgTheta(),
                         config.getBgLength(), config.getBgWidth(), 
                         config.getBgNx(), config.getBgNy(), config.getReynolds());
    DataStructure compMesh(config.getCompX0(), config.getCompY0(), config.getCompTheta(),
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
    
    cout << "Output directory: " << outputDir << endl;

    compMesh.MarkBoundaryDonor(bgMesh);
    bgMesh.HoleCutting(compMesh);
    bgMesh.MarkBoundaryPointType(BC_SPECIFIED);
    compMesh.MarkBoundaryPointType(INTERPOLATION_RECIEVER);
    
    bgMesh.StoreInternalInterpolatedPoints();
    compMesh.StoreInternalInterpolatedPoints();

    bgMesh.WriteTxtPointType(outputDir + "/bgPtType.txt");
    compMesh.WriteTxtPointType(outputDir + "/compPtType.txt");

    initialize(&bgMesh, &compMesh);
    initialize(&compMesh, &bgMesh);
    Solve(&bgMesh, &compMesh, config.getMaxIterations());
    GetOutput(&bgMesh, (outputDir + "/output_Upwind_bgMesh.dat").c_str());
    GetOutput(&compMesh, (outputDir + "/output_Upwind_compMesh.dat").c_str());
    // Also emit gnuplot data and a plotting script
    WriteGnuplotAll(bgMesh, compMesh, "NS contour for 2d plate", "lid_overset", outputDir);
    cout << "Solved overset mesh." << endl;

    return 0;
}
