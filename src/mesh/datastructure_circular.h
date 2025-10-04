#ifndef DATASTRUCTURE_CIRCULAR_H
#define DATASTRUCTURE_CIRCULAR_H

#include "mesh/datastructure.h"

// Derived class for circular/concentric mesh
class DataStructureCircular : public DataStructure {
   public:
    // Circular mesh specific parameters
    double xCenter, yCenter;         // Center of the circular domain
    double innerRadius, outerRadius; // Inner and outer radii for concentric circular hole
    int Nr, Ntheta;                  // Number of points in radial and angular directions
    double dr, dtheta;               // Radial and angular spacing

    // Constructor for circular mesh with concentric hole
    DataStructureCircular(double xC, double yC, double innerR, double outerR, 
                         size_t _Nr, size_t _Ntheta, double Reynolds = 100.0);
    virtual ~DataStructureCircular() = default;

    // Override virtual methods from base class
    void GenerateMesh() override;
    void StoreInternalInterpolatedPoints() override;
    size_t GetElementNumber(size_t i, size_t j) const override;
    size_t GetPointNumber(size_t i, size_t j) const override;
    std::tuple<size_t, size_t> GetijFromPointNumber(size_t iPoint) const override;

private:
    void GenerateCircularPointCoordinates();
    void GenerateCircularElements();
    void GenerateCircularBoundaries();
};

#endif // DATASTRUCTURE_CIRCULAR_H
