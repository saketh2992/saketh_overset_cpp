#ifndef SOLVER_H
#define SOLVER_H

#include "mesh/datastructure.h"

// Initialization functions
void CopyNewtoOld(DataStructure *rect);
void InterpolateCells(DataStructure *rect, DataStructure *oversetMesh, int k);
void ApplyBC(DataStructure *rect, int k, DataStructure *oversetMesh);
void LinearInterpolation(DataStructure *rect);
void initialize(DataStructure *rect, DataStructure *oversetMesh);

// Flux calculation functions
void SimpleUpwind(DataStructure *rect, double *Fc, double *ap_c, int i, int j, int k);
void Quick(DataStructure *rect, double *Fc, double *ap_c, int i, int j, int k);
void DiffusiveFlux(DataStructure *rect, double *Fd, double *ap_d, int i, int j, int k);
void UpdateFlux(DataStructure *rect);

// Solving functions
double SolveUV(DataStructure *rect, DataStructure *oversetMesh, int k);
double SolveP(DataStructure *rect, DataStructure *oversetMesh, int k, const double RELAX);
int CorrectVelocity(DataStructure *rect);

// Main solving functions
void ImplicitSolve(DataStructure *rect, DataStructure *oversetMesh, int outerIter, int& activePointsBg, int& activePointsComp);
bool ConvergenceCheck(DataStructure *rect, int count, int activePoints);
void Solve(DataStructure *rect, DataStructure *oversetMesh, int maxIterations = 60000);

// Utility functions
inline bool OuterIterLogFreq(const int iter);

#endif // SOLVER_H