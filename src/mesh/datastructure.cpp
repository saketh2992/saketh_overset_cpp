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

// Base class constructor - initializes common parameters
DataStructure::DataStructure() {
    Dim = 2;
    nVar = 3;

    /* Physical constants*/
    ulid = 1.0;
    rho = 1.0;

    Init.resize(nVar, 0.0);
    residual.resize(nVar);
}

// Common initialization for arrays that depend on numberOfPoints
void DataStructure::InitializeCommonArrays() {
    interpolationStencil.resize(numberOfPoints);
    interpolationCoeffs.resize(numberOfPoints);

    /* Point coordinates */
    pointCoordinates.resize(Dim);
    for (unsigned short iDim = 0; iDim < Dim; ++iDim) {
        pointCoordinates[iDim].resize(numberOfPoints);
    }
    pointType.resize(numberOfPoints, CALCULATED);
}

// Build element connectivity and bounding boxes - called by derived classes after GenerateMesh
void DataStructure::BuildElementConnectivity() {
    elementBBox.resize(numberOfElements);
    
    for (size_t iElement = 0; iElement < numberOfElements; ++iElement) {
        /* Element cartesian aligned Bounding Box */
        elementBBox[iElement].assign({su2double_highest, su2double_highest, su2double_lowest, su2double_lowest});
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

// Build neighbor relationships - called by derived classes
void DataStructure::BuildNeighborPointsOfPoint() {
    neighborPointsOfPoint.resize(numberOfPoints);
    
    for (size_t iElement = 0; iElement < numberOfElements; ++iElement) {
        size_t pointA, pointB;
        for (size_t iPoint = 0; iPoint < elementConnectivity[iElement].size(); iPoint++) {
            pointA = elementConnectivity[iElement][iPoint];
            pointB = elementConnectivity[iElement][(iPoint + 1) % 4];
            neighborPointsOfPoint[pointA].push_back(pointB);
            neighborPointsOfPoint[pointB].push_back(pointA);
        }
    }
    
    /* Remove duplicates from the neighboring point lists */
    vector<size_t>::iterator vecIt;
    for (size_t iPoint = 0; iPoint < numberOfPoints; ++iPoint) {
        sort(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());
        vecIt = unique(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());
        neighborPointsOfPoint[iPoint].resize(vecIt - neighborPointsOfPoint[iPoint].begin());
    }
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
