#include "mesh/datastructure_rectangular.h"
#include "util/utilities.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>

using namespace std;

DataStructureRectangular::DataStructureRectangular(double xOrigin, double yOrigin, double axisTheta, 
                                                 double lengthX, double lengthY, size_t _Nx, size_t _Ny, 
                                                 double Reynolds) : DataStructure() {
    // Set rectangular mesh specific parameters
    xO = xOrigin;
    yO = yOrigin;
    theta = axisTheta;
    Lx = lengthX;
    Ly = lengthY;
    Nx = _Nx;
    Ny = _Ny;
    dx = Lx / Nx;
    dy = Ly / Ny;

    numberOfPoints = (Nx + 2) * (Ny + 2);
    numberOfElements = (Nx + 1) * (Ny + 1);
    volp = dx * dy;

    /* Calculate physical constants*/
    nu = ulid * Lx / Reynolds;

    /* Courant Number < 1*/
    dt = 0.2 * min(dx, dy) / ulid;
    
    /* Mesh details*/
    cout << "Rectangular Mesh Details: " << xO << " " << yO << " " << Nx << " " << Ny 
         << " " << Lx << " " << Ly << " " << theta << " dT:" << dt << endl;

    // Initialize variable arrays
    Var.resize(nVar);
    VarOld.resize(nVar);
    for (int i = 0; i < nVar; i++) {
        Var[i].resize(Nx + 2);
        VarOld[i].resize(Nx + 2);
        for (int j = 0; j < (Nx + 2); j++) {
            Var[i][j].resize(Ny + 2, 0.0);
            VarOld[i][j].resize(Ny + 2, 0.0);
        }
    }

    Ff.resize(2 * Dim);
    for (int k = 0; k < (2 * Dim); k++) {
        Ff[k].resize(Nx + 2);
        for (int j = 0; j < (Nx + 2); j++) {
            Ff[k][j].resize(Ny + 2, 0.0);
        }
    }

    /* Normal vectors for faces. No of faces per cell = 4  i=0,1,2,3 -> E,N,W,S*/
    Sf.resize(2 * Dim);
    for (int k = 0; k < (2 * Dim); k++) {
        Sf[k].resize(Dim);
    }
    for (int k = 0; k < (2*Dim); ++k) {
        Sf[k][0] = cos(theta + double(k) * 90.0 * PI/180.0); // i comp
        Sf[k][1] = sin(theta + double(k) * 90.0 * PI/180.0); // j comp
    }

    // Generate mesh using virtual function
    GenerateMesh();
}

// Override: Generate rectangular mesh
void DataStructureRectangular::GenerateMesh() {
    InitializeCommonArrays();
    GenerateRectangularPointCoordinates();
    GenerateRectangularElements();
    BuildElementConnectivity();
    BuildNeighborPointsOfPoint();
    GenerateRectangularBoundaries();
    GenerateADT();
}

void DataStructureRectangular::GenerateRectangularPointCoordinates() {
    for (size_t j = 0; j < (Ny + 2); ++j) {
        for (size_t i = 0; i < (Nx + 2); i++) {
            size_t iPoint = GetPointNumber(i, j);
            /* cell centered so adding dx/2 and dy/2 to point coordinates*/
            pointCoordinates[0][iPoint] = dx / 2 + xO + (i * dx) * cos(theta) - (j * dy) * sin(theta);
            pointCoordinates[1][iPoint] = dy / 2 + yO + (i * dx) * sin(theta) + (j * dy) * cos(theta);
        }
    }
}

void DataStructureRectangular::GenerateRectangularElements() {
    elementConnectivity.resize(numberOfElements);
    for (size_t j = 0; j < Ny + 1; ++j) {
        for (size_t i = 0; i < Nx + 1; ++i) {
            size_t iElement = GetElementNumber(i, j);
            elementConnectivity[iElement].assign({
                j * (Nx + 2) + i, 
                j * (Nx + 2) + i + 1, 
                (j + 1) * (Nx + 2) + i + 1, 
                (j + 1) * (Nx + 2) + i
            });
        }
    }
}

void DataStructureRectangular::GenerateRectangularBoundaries() {
    boundaryPoints.resize(4); /* bottom, right, top, left*/
    for (size_t iBoundPoint = 0; iBoundPoint < (Nx + 2); iBoundPoint++) {
        boundaryPoints[0].push_back(iBoundPoint);                       /*bottom*/
        boundaryPoints[2].push_back((Ny + 1) * (Nx + 2) + iBoundPoint); /*top*/
    }
    for (size_t iBoundPoint = 0; iBoundPoint < (Ny + 2); iBoundPoint++) {
        boundaryPoints[3].push_back(iBoundPoint * (Nx + 2));            /*left*/
        boundaryPoints[1].push_back((Nx + 2) * iBoundPoint + (Nx + 1)); /*right*/
    }
    reverse(boundaryPoints[2].begin(), boundaryPoints[2].end());
    reverse(boundaryPoints[3].begin(), boundaryPoints[3].end());
}

// Override: Store internal interpolated points for rectangular mesh
void DataStructureRectangular::StoreInternalInterpolatedPoints() {
    for (int i = 1; i < Nx + 1; i++) {
        for (int j = 1; j < Ny + 1; j++) {
            size_t point = GetPointNumber(i, j);
            if (pointType[point] == INTERPOLATION_RECIEVER) {
                interpolatedPoints.push_back(point);
            }
        }
    }
}

// Override: Rectangular mesh indexing functions
size_t DataStructureRectangular::GetElementNumber(size_t i, size_t j) const { 
    return j * (Nx + 1) + i; 
}

size_t DataStructureRectangular::GetPointNumber(size_t i, size_t j) const { 
    return j * (Nx + 2) + i; 
}

std::tuple<size_t, size_t> DataStructureRectangular::GetijFromPointNumber(size_t iPoint) const { 
    return {iPoint % (Nx + 2), iPoint / (Nx + 2)}; 
}
