#include <iostream>
#include "constants.h"
#include "utilities.h"
#include "adt.h"
#include "datastructure.h"
#include "solver.h"
#include "output.h"

using namespace std;
int main() {
    const int mulFac = 5;
    DataStructure bgMesh(0.0, 0.0, 0.0 * PI / 180.0 , 1.0, 1.0, 10*mulFac, 10*mulFac);
    DataStructure compMesh(0.42, 0.42, 0.0 * PI / 180.0, 0.4, 0.4, 6*mulFac, 6*mulFac);

    compMesh.MarkBoundaryDonor(bgMesh);
    // bgMesh.HoleCutting(compMesh);
    bgMesh.MarkBoundaryPointType(BC_SPECIFIED);
    compMesh.MarkBoundaryPointType(INTERPOLATION_RECIEVER);
    
    bgMesh.StoreInternalInterpolatedPoints();
    compMesh.StoreInternalInterpolatedPoints();

    bgMesh.WriteTxtPointType("bgPtType.txt");
    compMesh.WriteTxtPointType("compPtType.txt");

    initialize(&bgMesh, &compMesh);
    initialize(&compMesh, &bgMesh);
    Solve(&bgMesh, &compMesh);
    GetOutput(&bgMesh, "output_Upwind_bgMesh.dat");
    GetOutput(&compMesh, "output_Upwind_compMesh.dat");
    cout << "Solved overset mesh." << endl;

    return 0;
}