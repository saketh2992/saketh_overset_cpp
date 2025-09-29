#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H

#include <vector>
#include <array>
#include <tuple>
#include <string>
#include "constants.h"
#include "adt.h"

class DataStructure {
   public:
    int Dim;
    int nVar;

    double xO, yO, theta;
    double dx, dy, Lx, Ly;
    int Nx, Ny, numberOfPoints, numberOfElements;
    double dt;
    double volp;

    double ulid;
    double nu;
    double rho;

    std::vector<std::vector<std::vector<double> > > Var;
    std::vector<std::vector<std::vector<double> > > VarOld;
    std::vector<std::vector<std::vector<double> > > Ff;
    std::vector<std::vector<double> > Sf;

    std::vector<double> Init;
    std::vector<double> residual;

    std::vector<unsigned short> pointType;
    std::vector<std::vector<double> > pointCoordinates;
    std::vector<std::vector<size_t> > neighborPointsOfPoint;
    std::vector<std::vector<size_t> > elementConnectivity;
    std::vector<std::vector<double> > elementBBox;

    std::vector<std::vector<size_t> > interpolationStencil;
    std::vector<std::vector<double> > interpolationCoeffs;

    std::vector<size_t> interpolatedPoints;

    /* bottom, right, top, left*/
    std::vector<std::vector<size_t> > boundaryPoints;
    /* bottom, right, top, left*/
    std::vector<double> boundaryCondition;

    ADT adtBoundingBox;

    // Constructor
    DataStructure(double xOrigin, double yOrigin, double axisTheta, double lengthX, double lengthY, size_t _Nx, size_t _Ny);
    virtual ~DataStructure() = default;

    // Inline utility functions
    size_t GetElementNumber(size_t i, size_t j) const;
    size_t GetPointNumber(size_t i, size_t j) const;
    std::tuple<size_t, size_t> GetijFromPointNumber(size_t iPoint) const;

    // Public methods
    void GenerateADT();
    std::array<su2double, SU2_BBOX_SIZE> GetPointBBox(const size_t pointNumber) const;
    std::vector<size_t> GetValidDonorBBox(const DataStructure &oversetMesh, const size_t iPoint);
    void StoreInternalInterpolatedPoints();
    void MarkBoundaryDonor(DataStructure &backgroundMesh);
    void HoleCutting(const DataStructure& oversetMesh);
    void MarkBoundaryPointType(const unsigned short ptType);
    void WriteTxtPointType(std::string filename = "pointType.txt") const;
};

#endif // DATASTRUCTURE_H