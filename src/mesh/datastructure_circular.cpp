#include "mesh/datastructure_circular.h"
#include "util/utilities.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>

using namespace std;

DataStructureCircular::DataStructureCircular(double xC, double yC, double innerR, double outerR, 
                                           size_t _Nr, size_t _Ntheta, double Reynolds) : DataStructure() {
    // Set circular mesh specific parameters
    xCenter = xC;
    yCenter = yC;
    innerRadius = innerR;
    outerRadius = outerR;
    Nr = _Nr;
    Ntheta = _Ntheta;
    
    dr = (outerRadius - innerRadius) / Nr;
    dtheta = 2.0 * PI / Ntheta;

    numberOfPoints = (Nr + 2) * Ntheta;
    numberOfElements = (Nr + 1) * Ntheta;
    volp = 0.5 * dr * (2.0 * innerRadius + dr) * dtheta;

    /* Calculate physical constants*/
    nu = ulid * (outerRadius - innerRadius) / Reynolds;
    dt = 0.2 * min(dr, innerRadius * dtheta) / ulid;
    
    cout << "Circular Mesh Details: Center(" << xCenter << ", " << yCenter << ") "
         << "Inner R: " << innerRadius << " Outer R: " << outerRadius 
         << " Nr: " << Nr << " Ntheta: " << Ntheta << " dT:" << dt << endl;

    // Initialize variable arrays
    Var.resize(nVar);
    VarOld.resize(nVar);
    for (int i = 0; i < nVar; i++) {
        Var[i].resize(Nr + 2);
        VarOld[i].resize(Nr + 2);
        for (int j = 0; j < (Nr + 2); j++) {
            Var[i][j].resize(Ntheta, 0.0);
            VarOld[i][j].resize(Ntheta, 0.0);
        }
    }
    
    Ff.resize(2 * Dim);
    for (int k = 0; k < (2 * Dim); k++) {
        Ff[k].resize(Nr + 2);
        for (int j = 0; j < (Nr + 2); j++) {
            Ff[k][j].resize(Ntheta, 0.0);
        }
    }

    /* Face normal vectors for circular geometry */
    Sf.resize(2 * Dim);
    for (int k = 0; k < (2 * Dim); k++) {
        Sf[k].resize(Dim);
    }

    // Generate mesh using virtual function
    GenerateMesh();
}

// Override: Generate circular mesh
void DataStructureCircular::GenerateMesh() {
    InitializeCommonArrays();
    GenerateCircularPointCoordinates();
    GenerateCircularElements();
    BuildElementConnectivity();
    BuildNeighborPointsOfPoint();
    GenerateCircularBoundaries();
    GenerateADT();
}

void DataStructureCircular::GenerateCircularPointCoordinates() {
    for (size_t j = 0; j < (Nr + 2); ++j) {  // radial direction
        for (size_t i = 0; i < Ntheta; ++i) {  // angular direction
            size_t iPoint = GetPointNumber(i, j);
            
            // Calculate radius for this radial layer
            // j=0: inner ghost, j=1 to Nr: actual points, j=Nr+1: outer ghost
            double r = innerRadius + (j - 0.5) * dr;  // cell-centered
            
            // Calculate angle for this angular position
            double theta = i * dtheta + dtheta / 2.0;  // cell-centered
            
            // Convert polar to Cartesian coordinates
            pointCoordinates[0][iPoint] = xCenter + r * cos(theta);
            pointCoordinates[1][iPoint] = yCenter + r * sin(theta);
        }
    }
}

void DataStructureCircular::GenerateCircularElements() {
    elementConnectivity.resize(numberOfElements);
    for (size_t j = 0; j < Nr + 1; ++j) {  // radial layers
        for (size_t i = 0; i < Ntheta; ++i) {  // angular divisions
            size_t iElement = GetElementNumber(i, j);
            
            // Quadrilateral element connectivity (counterclockwise)
            size_t i_next = (i + 1) % Ntheta;  // wrap around in angular direction
            
            elementConnectivity[iElement].assign({
                j * Ntheta + i,            // inner radius, current angle
                j * Ntheta + i_next,       // inner radius, next angle
                (j + 1) * Ntheta + i_next, // outer radius, next angle
                (j + 1) * Ntheta + i       // outer radius, current angle
            });
        }
    }
}

void DataStructureCircular::GenerateCircularBoundaries() {
    boundaryPoints.resize(2); /* inner boundary (0), outer boundary (1) */
    
    // Inner boundary (j = 0, all angular positions)
    for (size_t i = 0; i < Ntheta; ++i) {
        boundaryPoints[0].push_back(GetPointNumber(i, 0));
    }
    
    // Outer boundary (j = Nr+1, all angular positions)
    for (size_t i = 0; i < Ntheta; ++i) {
        boundaryPoints[1].push_back(GetPointNumber(i, Nr + 1));
    }
}

// Override: Store internal interpolated points for circular mesh
void DataStructureCircular::StoreInternalInterpolatedPoints() {
    for (int j = 1; j < Nr + 1; j++) {
        for (int i = 0; i < Ntheta; i++) {
            size_t point = GetPointNumber(i, j);
            if (pointType[point] == INTERPOLATION_RECIEVER) {
                interpolatedPoints.push_back(point);
            }
        }
    }
}

// Override: Circular mesh indexing functions
size_t DataStructureCircular::GetElementNumber(size_t i, size_t j) const { 
    return j * Ntheta + i; 
}

size_t DataStructureCircular::GetPointNumber(size_t i, size_t j) const { 
    return j * Ntheta + i; 
}

std::tuple<size_t, size_t> DataStructureCircular::GetijFromPointNumber(size_t iPoint) const { 
    return {iPoint % Ntheta, iPoint / Ntheta}; 
}
