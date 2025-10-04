#ifndef DATASTRUCTURE_RECTANGULAR_H
#define DATASTRUCTURE_RECTANGULAR_H

#include "mesh/datastructure.h"

// Derived class for rectangular mesh
class DataStructureRectangular : public DataStructure {
   public:
    // Rectangular mesh specific parameters
    double xO, yO, theta;
    double dx, dy, Lx, Ly;
    int Nx, Ny;

    // Constructor for rectangular mesh
    DataStructureRectangular(double xOrigin, double yOrigin, double axisTheta, 
                           double lengthX, double lengthY, size_t _Nx, size_t _Ny, 
                           double Reynolds = 100.0);
    virtual ~DataStructureRectangular() = default;

    // Override virtual methods from base class
    void GenerateMesh() override;
    void StoreInternalInterpolatedPoints() override;
    size_t GetElementNumber(size_t i, size_t j) const override;
    size_t GetPointNumber(size_t i, size_t j) const override;
    std::tuple<size_t, size_t> GetijFromPointNumber(size_t iPoint) const override;

private:
    void GenerateRectangularPointCoordinates();
    void GenerateRectangularElements();
    void GenerateRectangularBoundaries();
};

#endif // DATASTRUCTURE_RECTANGULAR_H
