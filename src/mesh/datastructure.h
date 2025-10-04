#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H

#include <vector>
#include <array>
#include <tuple>
#include <string>
#include "util/constants.h"
#include "mesh/adt.h"

// Base class for all mesh data structures
class DataStructure {
   public:
    int Dim;
    int nVar;

    // Common mesh parameters
    int numberOfPoints, numberOfElements;
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

    std::vector<std::vector<size_t> > boundaryPoints;
    std::vector<double> boundaryCondition;

    ADT adtBoundingBox;

    // Constructor
    DataStructure();
    virtual ~DataStructure() = default;

    // Virtual methods to be implemented by derived classes
    virtual void GenerateMesh() = 0;  // Pure virtual - must be implemented by derived classes
    virtual size_t GetElementNumber(size_t i, size_t j) const = 0;
    virtual size_t GetPointNumber(size_t i, size_t j) const = 0;
    virtual std::tuple<size_t, size_t> GetijFromPointNumber(size_t iPoint) const = 0;

    // Common public methods
    void GenerateADT();
    std::array<su2double, SU2_BBOX_SIZE> GetPointBBox(const size_t pointNumber) const;
    std::vector<size_t> GetValidDonorBBox(const DataStructure &oversetMesh, const size_t iPoint);
    virtual void StoreInternalInterpolatedPoints() = 0;  // Virtual - mesh-specific implementation
    void MarkBoundaryDonor(DataStructure &backgroundMesh);
    void HoleCutting(const DataStructure& oversetMesh);
    void MarkBoundaryPointType(const unsigned short ptType);
    void WriteTxtPointType(std::string filename = "pointType.txt") const;

protected:
    // Protected helper methods for common initialization
    void InitializeCommonArrays();
    void BuildElementConnectivity();
    void BuildNeighborPointsOfPoint();
};

#endif // DATASTRUCTURE_H