#include "mesh/datastructure.h"
#include "util/utilities.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>

using namespace std;

DataStructure::DataStructure(double xOrigin, double yOrigin, double axisTheta, double lengthX, double lengthY, size_t _Nx, size_t _Ny, double Reynolds) {
    Dim = 2;
    nVar = 3;

    xO = xOrigin, yO = yOrigin, theta = axisTheta;
    Lx = lengthX, Ly = lengthY;
    Nx = _Nx, Ny = _Ny;
    dx = Lx / Nx;
    dy = Ly / Ny;

    numberOfPoints = (Nx + 2) * (Ny + 2);
    numberOfElements = (Nx + 1) * (Ny + 1);
    volp = dx * dy;

    /* Physical constants*/
    ulid = 1.0;
    rho = 1.0;
    nu = ulid * Lx / Reynolds;  // Calculate kinematic viscosity from Reynolds number

    /* Courant Number < 1*/
    dt = 0.2 * min(dx, dy)/ ulid; // taking 20% of the limit from Co < 1
    // dt = 0.001;
    
    /* Mesh details*/
    cout << "Mesh Details: " << xO << " " << yO << " " << Nx << " " << Ny << " " << Lx << " " << Ly << " " << theta << " dT:" << dt << endl;

    Init.resize(nVar, 0.0);
    Var.resize(nVar);
    VarOld.resize(nVar);

    /* normal vector for faces. No of faces per cell = 4  i=0,1,2,3 -> E,N,W,S*/
    Sf.resize(2 * Dim);
    for (int k = 0; k < (2 * Dim); k++) {
        Sf[k].resize(Dim);
    }
    /* east */
    for (int k = 0; k < (2*Dim); ++k) {
        Sf[k][0] = cos(theta + double(k) * 90.0 * PI/180.0); // i comp
        Sf[k][1] = sin(theta + double(k) * 90.0 * PI/180.0); // j comp
    }

    interpolationStencil.resize(numberOfPoints);
    interpolationCoeffs.resize(numberOfPoints);

    for (int i = 0; i < nVar; i++) {
        Var[i].resize(Nx + 2);
        VarOld[i].resize(Nx + 2);
        for (int j = 0; j < (Nx + 2); j++) {
            Var[i][j].resize(Ny + 2, 0.0);
            VarOld[i][j].resize(Ny + 2, 0.0);
        }
    }
    residual.resize(nVar);
    Ff.resize(2 * Dim);  // No of faces per cell = 4  i=0,1,2,3 -> E,N,W,S
    for (int k = 0; k < (2 * Dim); k++) {
        Ff[k].resize(Nx + 2);
        for (int j = 0; j < (Nx + 2); j++) {
            Ff[k][j].resize(Ny + 2, 0.0);
        }
    }

    /* requirements for overset mesh */
    /* Point coordinates */
    pointCoordinates.resize(Dim);
    for (unsigned short iDim = 0; iDim < Dim; ++iDim) {
        pointCoordinates[iDim].resize(numberOfPoints);
    }
    pointType.resize(numberOfPoints, CALCULATED);
    size_t iPoint = numberOfPoints;
    for (size_t j = 0; j < (Ny + 2); ++j) {
        for (size_t i = 0; i < (Nx + 2); i++) {
            iPoint = GetPointNumber(i, j);
            /* cell centered so adding dx/2 and dy/2 to point coordinates*/
            pointCoordinates[0][iPoint] = dx / 2 + xO + (i * dx) * cos(theta) - (j * dy) * sin(theta); /* global x axis*/
            pointCoordinates[1][iPoint] = dy / 2 + yO + (i * dx) * sin(theta) + (j * dy) * cos(theta); /* global y axis*/
        }
    }

    /* Element connectivity and it's bounding box */
    elementConnectivity.resize(numberOfElements);
    elementBBox.resize(numberOfElements);
    size_t iElement = numberOfElements;
    for (size_t j = 0; j < Ny + 1; ++j) {
        for (size_t i = 0; i < Nx + 1; ++i) {
            iElement = GetElementNumber(i, j);
            elementConnectivity[iElement].assign({j * (Nx + 2) + i, j * (Nx + 2) + i + 1, (j + 1) * (Nx + 2) + i + 1, (j + 1) * (Nx + 2) + i});

            /* Element cartesian aligned Bounding Box */
            elementBBox[iElement].assign({su2double_highest, su2double_highest, su2double_lowest, su2double_lowest}); /* {Xmin, Ymin, Xmax, Ymax}*/
            for (size_t iPoint = 0; iPoint < elementConnectivity[iElement].size(); ++iPoint) {
                double x = pointCoordinates[0][elementConnectivity[iElement][iPoint]];
                double y = pointCoordinates[1][elementConnectivity[iElement][iPoint]];
                elementBBox[iElement][0] = min(elementBBox[iElement][0], x);
                elementBBox[iElement][1] = min(elementBBox[iElement][1], y);
                elementBBox[iElement][2] = max(elementBBox[iElement][2], x);
                elementBBox[iElement][3] = max(elementBBox[iElement][3], y);
            }
        }
    }

    /* Neighbour points of point */
    neighborPointsOfPoint.resize(numberOfPoints);
    for (size_t j = 0; j < Ny + 1; ++j) {
        for (size_t i = 0; i < Nx + 1; ++i) {
            iElement = GetElementNumber(i, j);
            size_t pointA, pointB;
            for (size_t iPoint = 0; iPoint < elementConnectivity[iElement].size(); iPoint++) {
                pointA = elementConnectivity[iElement][iPoint];
                pointB = elementConnectivity[iElement][(iPoint + 1) % 4];
                neighborPointsOfPoint[pointA].push_back(pointB);
                neighborPointsOfPoint[pointB].push_back(pointA);
            }
        }
    }
    /* remove duplicates from the neighboring point lists*/
    iPoint = numberOfPoints;
    vector<size_t>::iterator vecIt;
    for (size_t j = 0; j < Ny + 2; ++j) {
        for (size_t i = 0; i < Nx + 2; ++i) {
            iPoint = GetPointNumber(i, j);
            /* sort neighboring points for each point */
            sort(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());

            /* uniquify list of neighboring points */
            vecIt = unique(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());

            /* adjust size of vector */
            neighborPointsOfPoint[iPoint].resize(vecIt - neighborPointsOfPoint[iPoint].begin());
        }
    }

    /* Boundary Surfaces*/
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

    GenerateADT();
}

size_t DataStructure::GetElementNumber(size_t i, size_t j) const { 
    return j * (Nx + 1) + i; 
}

size_t DataStructure::GetPointNumber(size_t i, size_t j) const { 
    return j * (Nx + 2) + i; 
}

std::tuple<size_t, size_t> DataStructure::GetijFromPointNumber(size_t iPoint) const { 
    return {iPoint % (Nx + 2), iPoint / (Nx + 2)}; 
}

void DataStructure::GenerateADT() {
    array<su2double, SU2_BBOX_SIZE> bboxCoords{};

    /* modifying node insertion to balance the tree */
    vector<size_t> elementOrderADT;
    elementOrderADT.resize(numberOfElements);
    iota(elementOrderADT.begin(), elementOrderADT.end(), 0);
    // shuffle(elementOrderADT.begin(), elementOrderADT.end(), std::mt19937{std::random_device{}()});
    shuffle(elementOrderADT.begin(), elementOrderADT.end(), std::mt19937{12337});

    for (size_t randomIndex = 0; randomIndex < numberOfElements; ++randomIndex) {
        size_t LocalIndex = elementOrderADT[randomIndex];
        for (unsigned short i = 0; i < SU2_BBOX_SIZE; ++i) {
            bboxCoords[i] = elementBBox[LocalIndex][i];
        }
        node_adt *tempNode = new node_adt(LocalIndex, bboxCoords);
        adtBoundingBox.insertNode(tempNode);
    }
    cout << "Max heirarchy : " << adtBoundingBox.treeHeirarchy << endl;
}

array<su2double, SU2_BBOX_SIZE> DataStructure::GetPointBBox(const size_t pointNumber) const {
    /* small value to generate a bounding box centered around a point */
    su2double epsilon = std::numeric_limits<su2double>::epsilon();
    array<su2double, 3> Coords{0.0, 0.0, 0.0};
    array<su2double, SU2_BBOX_SIZE> pointBBoxCoords{};
    for (unsigned short iDim = 0; iDim < Dim; ++iDim) {
        Coords[iDim] = pointCoordinates[iDim][pointNumber];
        pointBBoxCoords[iDim] = Coords[iDim] - epsilon;
        pointBBoxCoords[iDim + 2] = Coords[iDim] + epsilon;
    }
    return pointBBoxCoords;
}

vector<size_t> DataStructure::GetValidDonorBBox(const DataStructure &oversetMesh,const size_t iPoint) {
    const double TOLERANCE = 1e-6;
    const auto pointBBox = GetPointBBox(iPoint);
    const auto intersectingBBoxList = oversetMesh.adtBoundingBox.searchADT(pointBBox);
    
    bool validInterpolation = true;
    vector<size_t> validBBox;
    vector<double> tempInterpolationCoeffs;
    for (auto testBBox : intersectingBBoxList) {
        auto donorStencil = oversetMesh.elementConnectivity[testBBox];
        vector<vector<double> > donorCoords(Dim);
        vector<size_t> tempInterpolationStencil;
        for (size_t iDonor = 0; iDonor < donorStencil.size(); ++iDonor) {
            tempInterpolationStencil.push_back(donorStencil[iDonor]);
            for (unsigned short iDim = 0; iDim < Dim; ++iDim) {
                donorCoords[iDim].push_back(oversetMesh.pointCoordinates[iDim][donorStencil[iDonor]]);
            }
        }
        vector<double> receiverCoords;
        for (unsigned short iDim = 0; iDim < Dim; ++iDim) {
            receiverCoords.push_back(pointCoordinates[iDim][iPoint]);
        }
        tempInterpolationCoeffs = interpolation_coeff(donorCoords, receiverCoords);
        
        validInterpolation = true; // all coeff should lie between 0 and 1
        for (auto coeff : tempInterpolationCoeffs) {
            validInterpolation = validInterpolation && ((abs(coeff-1.0) <= TOLERANCE) || (coeff < 1.0)) && ((abs(coeff-0.0) <= TOLERANCE) || (coeff > 0.0));
        }
        if (validInterpolation == true) {
            interpolationStencil[iPoint] = tempInterpolationStencil;
            interpolationCoeffs[iPoint] = tempInterpolationCoeffs;
            validBBox.push_back(testBBox);
            break;
        }
        else {
            validInterpolation = true;
        }
    }
    if (validBBox.size() == 0) {
        // cerr << "Point: " << iPoint << " No interpolation stencil found. What to do ??" << endl;
    }
    if (validBBox.size() > 1) {
        cerr << "Multiple valid donor found. NOT POSSIBLE." << endl;
    }
    return validBBox;
}

void DataStructure::StoreInternalInterpolatedPoints () {
    for (int i = 1; i < Nx + 1; i++) {
        for (int j = 1; j < Ny + 1; j++) {
            size_t point = GetPointNumber(i, j);
            if (pointType[point] == INTERPOLATION_RECIEVER) {
                interpolatedPoints.push_back(point);
            }
        }
    }
}

void DataStructure::MarkBoundaryDonor(DataStructure &backgroundMesh) {
    /* iterate over all the boundaries of current overset mesh to mark donors on background mesh */
    for (size_t iBound = 0; iBound < boundaryPoints.size(); ++iBound) {
        /* iterate over all points in iBound*/
        for (size_t iBoundPoint = 0; iBoundPoint < boundaryPoints[iBound].size(); ++iBoundPoint) {
            size_t point = boundaryPoints[iBound][iBoundPoint];
            pointType[point] = INTERPOLATION_RECIEVER;
            const auto BBox = GetValidDonorBBox(backgroundMesh, point);
            if (BBox.size() == 0) {
                cerr << "Point: " << point << " No interpolation stencil found for boundary of overset mesh." << endl;
            }
            for (auto iPoint : backgroundMesh.elementConnectivity[BBox[0]]) {
                // cout << "Marking iPoint = " << iPoint << " as donor " << endl;
                backgroundMesh.pointType[iPoint] = INTERPOLATION_DONOR;
                /* mark all neighbors as required (for complete stencil of donor)*/
                for (auto neighborPoint : backgroundMesh.neighborPointsOfPoint[iPoint]) {
                    if (backgroundMesh.pointType[neighborPoint] == INTERPOLATION_DONOR) {
                        continue;
                    }
                    backgroundMesh.pointType[neighborPoint] = DONOR_BUFFER; /*donor buffer*/
                }
            }
        }
    }
}

void DataStructure::HoleCutting(const DataStructure& oversetMesh) {
    /* does a type of holecutting / removal of coarse cells where overlapping finer cells are present*/
    /*check removability of each point -> if neighbours can be interpolation cells and cell is bigger than overlapping cell*/
    for (size_t iPoint = 0; iPoint < numberOfPoints; ++iPoint) {
        /* skip interpolation donors and it's buffer */
        if (pointType[iPoint] == INTERPOLATION_DONOR || pointType[iPoint] == DONOR_BUFFER) {
            continue;
        }

        /* TODO: check if donors/overlapping cells are finer/smaller than iPoint. Assumed for now.*/
        bool neighborsCanBeInterpolated = true;
        for (auto neighborPoint : neighborPointsOfPoint[iPoint]) {
            if (pointType[neighborPoint] == INTERPOLATION_DONOR || pointType[neighborPoint] == DONOR_BUFFER) {
                neighborsCanBeInterpolated = false;
                break;
            }
            const auto BBox = GetValidDonorBBox(oversetMesh, neighborPoint);
            if (BBox.size() == 0) {
                neighborsCanBeInterpolated = false;
            }
        }
        if (neighborsCanBeInterpolated) {
            pointType[iPoint] = UNUSED;  // iPoint becomes unused
            for (auto neighborPoint : neighborPointsOfPoint[iPoint]) {
                if (pointType[neighborPoint] != CALCULATED) {
                    continue;
                }
                pointType[neighborPoint] = INTERPOLATION_RECIEVER;  // neighbors are marked as interpolated
            }
        }
    }
}

/* pType : 3 = interpolated, 5 = boundary condition specified */
void DataStructure::MarkBoundaryPointType(const unsigned short ptType) {
    /* marks all the boundaries (for component meshes marked as intepolated. Assumed component meshes dont overlap each other nor do they intersect boundaries) */
    for (size_t iBound = 0; iBound < boundaryPoints.size(); ++iBound) {
        /* iterate over all points in iBound*/
        for (size_t iBoundPoint = 0; iBoundPoint < boundaryPoints[iBound].size(); ++iBoundPoint) {
            size_t point = boundaryPoints[iBound][iBoundPoint];
            pointType[point] = ptType;
        }
    }
}

void DataStructure::WriteTxtPointType(string filename) const {
    ofstream pointTypeFile(filename);
    if (!pointTypeFile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
    } else {
        for (size_t iPoint = 0; iPoint < numberOfPoints; iPoint++) {
            // if (iPoint % (Nx) == 0 && iPoint !=0) {pointTypeFile << endl;}
            pointTypeFile << pointType[iPoint] << " ";
        }
        pointTypeFile.close();
    }
}
