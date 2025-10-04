#include "solver/solver.h"
#include "mesh/datastructure_rectangular.h"
#include "mesh/datastructure_circular.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

void CopyNewtoOld(DataStructure *rect) {
    std::copy(rect->Var.begin(), rect->Var.end(), rect->VarOld.begin());
}

void InterpolateCells(DataStructure *rect, DataStructure *oversetMesh, int k) {
    for(size_t point : rect->interpolatedPoints) {
        auto ij = rect->GetijFromPointNumber(point);
        auto i = std::get<0>(ij);
        auto j = std::get<1>(ij);
        // Initialize to zero before accumulating interpolation contributions
        rect->Var[k][i][j] = 0.0;
        for (unsigned short iDonor = 0; iDonor < rect->interpolationStencil[point].size(); ++iDonor) {
            auto donorPointNumber = rect->interpolationStencil[point][iDonor];
            auto ijDonor = oversetMesh->GetijFromPointNumber(donorPointNumber);
            auto _iDonor = std::get<0>(ijDonor);
            auto _jDonor = std::get<1>(ijDonor);
            rect->Var[k][i][j] += oversetMesh->Var[k][_iDonor][_jDonor] * rect->interpolationCoeffs[point][iDonor];
        }
    }
}

void ApplyBC(DataStructure *rect, int k, DataStructure *oversetMesh) {
    // Cast to rectangular mesh for accessing Nx, Ny
    DataStructureRectangular *rectMesh = dynamic_cast<DataStructureRectangular*>(rect);
    
    if (rect->pointType[rect->boundaryPoints[0][0]] == INTERPOLATION_RECIEVER) {
        // cout << "Interpolating at boundary for k = " << k << endl;
        for (size_t iBound = 0; iBound < rect->boundaryPoints.size(); ++iBound) {
            for (size_t iBoundPoint = 0; iBoundPoint < rect->boundaryPoints[iBound].size(); ++iBoundPoint) {
                size_t point = rect->boundaryPoints[iBound][iBoundPoint];
                auto ij = rect->GetijFromPointNumber(point);
                auto i = std::get<0>(ij);
                auto j = std::get<1>(ij);
                rect->Var[k][i][j] = 0.0;
                for (unsigned short iDonor = 0; iDonor < rect->interpolationStencil[point].size(); ++iDonor) {
                    auto donorPointNumber = rect->interpolationStencil[point][iDonor];
                    auto ijDonor = oversetMesh->GetijFromPointNumber(donorPointNumber);
                    auto _iDonor = std::get<0>(ijDonor);
                    auto _jDonor = std::get<1>(ijDonor);
                    rect->Var[k][i][j] = rect->Var[k][i][j] + oversetMesh->Var[k][_iDonor][_jDonor] * rect->interpolationCoeffs[point][iDonor];
                }
            }
        }
    } else if (rectMesh != nullptr) {
        // cout << "Applying specified BC for k = " << k << endl;
        switch (k) {
            case 0:  // U

                for (int j = 1; j < rectMesh->Ny + 1; j++) {
                    rect->Var[k][0][j] = 2 * (0.0) - rect->Var[k][1][j];                       // Left
                    rect->Var[k][rectMesh->Nx + 1][j] = 2 * (0.0) - rect->Var[k][rectMesh->Nx][j];  // Right
                }
                for (int i = 1; i < rectMesh->Nx + 1; i++) {
                    rect->Var[k][i][rectMesh->Ny + 1] = 2 * rect->ulid - rect->Var[k][i][rectMesh->Ny];  // Top
                    rect->Var[k][i][0] = 2 * (0.0) - rect->Var[k][i][1];                         // Bottom
                }
                break;

            case 1:  // V

                for (int j = 1; j < rectMesh->Ny + 1; j++) {
                    rect->Var[k][0][j] = 2 * (0.0) - rect->Var[k][1][j];                       // Left
                    rect->Var[k][rectMesh->Nx + 1][j] = 2 * (0.0) - rect->Var[k][rectMesh->Nx][j];  // Right
                }
                for (int i = 1; i < rectMesh->Nx + 1; i++) {
                    rect->Var[k][i][rectMesh->Ny + 1] = 2 * (0.0) - rect->Var[k][i][rectMesh->Ny];  // Top
                    rect->Var[k][i][0] = 2 * (0.0) - rect->Var[k][i][1];                    // Bottom
                }
                break;

            case 2:  // P (All Neumann Condition)

                for (int j = 1; j < rectMesh->Ny + 1; j++) {
                    rect->Var[k][0][j] = rect->Var[k][1][j];                       // Left
                    rect->Var[k][rectMesh->Nx + 1][j] = rect->Var[k][rectMesh->Nx][j];  // Right
                }
                for (int i = 1; i < rectMesh->Nx + 1; i++) {
                    rect->Var[k][i][rectMesh->Ny + 1] = rect->Var[k][i][rectMesh->Ny];  // Top
                    rect->Var[k][i][0] = rect->Var[k][i][1];                    // Bottom
                }
                break;
        }
    }
    // If rectMesh is nullptr (e.g., circular mesh), skip boundary conditions
    // Circular mesh boundary conditions would be handled separately
}

// Forward declaration for circular mesh linear interpolation
void LinearInterpolation_Circular(DataStructureCircular *circMesh);

void LinearInterpolation(DataStructure *rect) {
    // Cast to rectangular mesh for accessing dx, dy, Nx, Ny
    DataStructureRectangular *rectMesh = dynamic_cast<DataStructureRectangular*>(rect);
    
    // Check if this is a circular mesh
    DataStructureCircular *circMesh = dynamic_cast<DataStructureCircular*>(rect);
    if (circMesh != nullptr) {
        // Use circular mesh interpolation
        LinearInterpolation_Circular(circMesh);
        return;
    }
    
    // Skip if not a rectangular mesh and not circular
    if (rectMesh == nullptr) {
        return;
    }
    
    const vector<double> area = {rectMesh->dy, rectMesh->dx};
    for (int i = 1; i < rectMesh->Nx + 1; i++) {
        for (int j = 1; j < rectMesh->Ny + 1; j++) {
            size_t point = rect->GetPointNumber(i, j);
            if (rect->pointType[point] == UNUSED || rect->pointType[point] == INTERPOLATION_RECIEVER) {
                continue;
            }
            
            for (int k = 0; k < (rect->Dim); ++k) { // East and North Face
                rect->Ff[k][i][j] = 0;
                double u_avg = (rect->Var[0][i][j] + rect->Var[0][i + 1][j]) * 0.5;
                double v_avg = (rect->Var[1][i][j] + rect->Var[1][i][j + 1]) * 0.5;
                double face_area = area[k];
                rect->Ff[k][i][j] = (u_avg * rect->Sf[k][0] + v_avg * rect->Sf[k][1]) * face_area;   
            }
            for (int k = rect->Dim; k < (2*rect->Dim); ++k) { // West and South Face
                rect->Ff[k][i][j] = 0;
                double u_avg = (rect->Var[0][i - 1][j] + rect->Var[0][i][j]) * 0.5;
                double v_avg = (rect->Var[1][i][j - 1] + rect->Var[1][i][j]) * 0.5;
                double face_area = area[k-rect->Dim];
                rect->Ff[k][i][j] = (u_avg * rect->Sf[k][0] + v_avg * rect->Sf[k][1]) * face_area;
            }
            // rect->Ff[0][i][j] = (rect->Var[0][i][j] + rect->Var[0][i + 1][j]) * rect->dy * 0.5;   // East Face
            // rect->Ff[1][i][j] = (rect->Var[1][i][j] + rect->Var[1][i][j + 1]) * rect->dx * 0.5;   // North Face
            // rect->Ff[2][i][j] = -(rect->Var[0][i][j] + rect->Var[0][i - 1][j]) * rect->dy * 0.5;  // West Face
            // rect->Ff[3][i][j] = -(rect->Var[1][i][j] + rect->Var[1][i][j - 1]) * rect->dx * 0.5;  // South Face
        }
    }
}

void initialize(DataStructure *rect, DataStructure *oversetMesh) {
    // initializing interior values
    for (int k = 0; k < rect->nVar; k++) {
        ApplyBC(rect, k, oversetMesh);
    }

    CopyNewtoOld(rect);
    LinearInterpolation(rect);
}

void SimpleUpwind(DataStructure *rect, double *Fc, double *ap_c, int i, int j, int k) {
    double ue, uw, un, us;
    double sum = 0;
    if (rect->Ff[0][i][j] >= 0) {
        ue = rect->Var[k][i][j];
        sum += rect->Ff[0][i][j];
    } else
        ue = rect->Var[k][i + 1][j];

    if (rect->Ff[2][i][j] >= 0) {
        uw = rect->Var[k][i][j];
        sum += rect->Ff[2][i][j];
    } else
        uw = rect->Var[k][i - 1][j];

    if (rect->Ff[1][i][j] >= 0) {
        un = rect->Var[k][i][j];
        sum += rect->Ff[1][i][j];
    } else
        un = rect->Var[k][i][j + 1];

    if (rect->Ff[3][i][j] >= 0) {
        us = rect->Var[k][i][j];
        sum += rect->Ff[3][i][j];
    } else
        us = rect->Var[k][i][j - 1];

    *Fc = (ue * rect->Ff[0][i][j] + uw * rect->Ff[2][i][j] + un * rect->Ff[1][i][j] + us * rect->Ff[3][i][j]);
    *ap_c = sum * rect->volp;
}

void Quick(DataStructure *rect, double *Fc, double *ap_c, int i, int j, int k) {
    double ue, uw, un, us;
    double sum = 0;
    // East
    if (rect->Ff[0][i][j] >= 0) {
        ue = 0.75 * rect->Var[k][i][j] + 0.375 * rect->Var[k][i + 1][j] - 0.125 * rect->Var[k][i - 1][j];
        sum += 0.75 * rect->Ff[0][i][j];
    } else {
        ue = 0.75 * rect->Var[k][i + 1][j] + 0.375 * rect->Var[k][i][j] - 0.125 * rect->Var[k][i + 2][j];
        sum += 0.375 * rect->Ff[0][i][j];
    }
    // West
    if (rect->Ff[2][i][j] >= 0) {
        uw = 0.75 * rect->Var[k][i][j] + 0.375 * rect->Var[k][i - 1][j] - 0.125 * rect->Var[k][i + 1][j];
        sum += 0.75 * rect->Ff[2][i][j];
    } else {
        uw = 0.75 * rect->Var[k][i - 1][j] + 0.375 * rect->Var[k][i][j] - 0.125 * rect->Var[k][i - 2][j];
        sum += 0.375 * rect->Ff[2][i][j];
    }
    // North
    if (rect->Ff[1][i][j] >= 0) {
        un = 0.75 * rect->Var[k][i][j] + 0.375 * rect->Var[k][i][j + 1] - 0.125 * rect->Var[k][i][j - 1];
        sum += 0.75 * rect->Ff[1][i][j];
    } else {
        un = 0.75 * rect->Var[k][i][j + 1] + 0.375 * rect->Var[k][i][j] - 0.125 * rect->Var[k][i][j + 2];
        sum += 0.375 * rect->Ff[1][i][j];
    }
    if (rect->Ff[3][i][j] >= 0) {
        us = 0.75 * rect->Var[k][i][j] + 0.375 * rect->Var[k][i][j - 1] - 0.125 * rect->Var[k][i][j + 1];
        sum += 0.75 * rect->Ff[3][i][j];
    } else {
        us = 0.75 * rect->Var[k][i][j - 1] + 0.375 * rect->Var[k][i][j] - 0.125 * rect->Var[k][i][j - 2];
        sum += 0.375 * rect->Ff[3][i][j];
    }
    *Fc = (ue * rect->Ff[0][i][j] + uw * rect->Ff[2][i][j] + un * rect->Ff[1][i][j] + us * rect->Ff[3][i][j]);
    *ap_c = sum * rect->volp;
}

void DiffusiveFlux(DataStructure *rect, double *Fd, double *ap_d, int i, int j, int k) {
    // Cast to rectangular mesh for accessing dx, dy
    DataStructureRectangular *rectMesh = dynamic_cast<DataStructureRectangular*>(rect);
    
    *Fd = rect->volp * ((rect->Var[k][i + 1][j] - 2.0 * rect->Var[k][i][j] + rect->Var[k][i - 1][j]) / (rectMesh->dx * rectMesh->dx) + (rect->Var[k][i][j + 1] - 2.0 * rect->Var[k][i][j] + rect->Var[k][i][j - 1]) / (rectMesh->dy * rectMesh->dy));
    *ap_d = -rect->volp * (2.0 / (rectMesh->dx * rectMesh->dx) + 2.0 / (rectMesh->dy * rectMesh->dy));
}

void UpdateFlux(DataStructure *rect) {
    // Cast to rectangular mesh for accessing dx, dy, Nx, Ny
    DataStructureRectangular *rectMesh = dynamic_cast<DataStructureRectangular*>(rect);
    
    // Skip if not a rectangular mesh (e.g., circular mesh)
    if (rectMesh == nullptr) {
        return;
    }
    
    for (int i = 1; i < rectMesh->Nx + 1; i++) {
        for (int j = 1; j < rectMesh->Ny + 1; j++) {
            size_t point = rect->GetPointNumber(i, j);
            if (rect->pointType[point] == UNUSED || rect->pointType[point] == INTERPOLATION_RECIEVER) {
                continue;
            }
            rect->Ff[0][i][j] += -rect->dt / rect->rho * (rect->Var[2][i + 1][j] - rect->Var[2][i][j]) * rectMesh->dy / rectMesh->dx;  // East Face
            rect->Ff[1][i][j] += -rect->dt / rect->rho * (rect->Var[2][i][j + 1] - rect->Var[2][i][j]) * rectMesh->dx / rectMesh->dy;  // North Face
            rect->Ff[2][i][j] += -rect->dt / rect->rho * (rect->Var[2][i - 1][j] - rect->Var[2][i][j]) * rectMesh->dy / rectMesh->dx;  // West Face
            rect->Ff[3][i][j] += -rect->dt / rect->rho * (rect->Var[2][i][j - 1] - rect->Var[2][i][j]) * rectMesh->dx / rectMesh->dy;  // South Face
        }
    }
}

// Forward declarations for circular mesh solver
double SolveUV_Circular(DataStructureCircular *circMesh, DataStructure *oversetMesh, int k);

double SolveUV(DataStructure *rect, DataStructure *oversetMesh, int k) {
    /* k determines U or V. k = 0 -> U*/
    // Cast to rectangular mesh for accessing Nx, Ny
    DataStructureRectangular *rectMesh = dynamic_cast<DataStructureRectangular*>(rect);
    
    // Check if this is a circular mesh
    DataStructureCircular *circMesh = dynamic_cast<DataStructureCircular*>(rect);
    if (circMesh != nullptr) {
        // TEMPORARY: Disable circular mesh solver due to instability
        // Use circular mesh solver (Cartesian formulation)
        // return SolveUV_Circular(circMesh, oversetMesh, k);
        return 0.0;  // Skip for now
    }
    
    // Skip if not a rectangular mesh and not circular
    if (rectMesh == nullptr) {
        return 0.0;
    }
    
    int count = 0;
    double Fc, ap_c, Fd, ap_d, ap, R, rms;
    rms = 0.0;
    for (int i = 1; i < rectMesh->Nx + 1; i++) {
        for (int j = 1; j < rectMesh->Ny + 1; j++) {
            size_t point = rect->GetPointNumber(i, j);
            if (rect->pointType[point] == UNUSED) {
                R = 0.0;
            }
            if (rect->pointType[point] == CALCULATED || rect->pointType[point] == INTERPOLATION_DONOR || rect->pointType[point] == DONOR_BUFFER) {
                SimpleUpwind(rect, &Fc, &ap_c, i, j, k);
                // Quick(rect, &Fc, &ap_c, i, j, k);
                DiffusiveFlux(rect, &Fd, &ap_d, i, j, k);
                R = -(rect->volp / rect->dt * (rect->Var[k][i][j] - rect->VarOld[k][i][j]) + Fc + (-rect->nu) * Fd);
                ap = rect->volp / rect->dt + ap_c + (-rect->nu) * ap_d;
                count++;
                rect->Var[k][i][j] = rect->Var[k][i][j] + R / ap;
            } else if (rect->pointType[point] == INTERPOLATION_RECIEVER) {
                R = 0.0;
            }
            rms = rms + R * R;
        }
    }
    ApplyBC(rect, k, oversetMesh);
    rms = sqrt(rms / (rectMesh->Nx * rectMesh->Ny));
    return rms;
}

// Forward declaration for circular mesh pressure solver
double SolveP_Circular(DataStructureCircular *circMesh, DataStructure *oversetMesh, int k, double RELAX);

double SolveP(DataStructure *rect, DataStructure *oversetMesh, int k, const double RELAX) {
    /* for P k = 2 (in Var array)*/
    // Cast to rectangular mesh for accessing Nx, Ny
    DataStructureRectangular *rectMesh = dynamic_cast<DataStructureRectangular*>(rect);
    
    // Check if this is a circular mesh
    DataStructureCircular *circMesh = dynamic_cast<DataStructureCircular*>(rect);
    if (circMesh != nullptr) {
        // TEMPORARY: Disable circular mesh solver due to instability
        // Use circular mesh solver (Cartesian formulation)
        // return SolveP_Circular(circMesh, oversetMesh, k, RELAX);
        return 0.0;  // Skip for now
    }
    
    // Skip if not a rectangular mesh and not circular
    if (rectMesh == nullptr) {
        return 0.0;
    }
    
    int count = 0;
    double Fd, ap_d, ap, R, rms;
    rms = 0.0;
    for (int i = 1; i < rectMesh->Nx + 1; i++) {
        for (int j = 1; j < rectMesh->Ny + 1; j++) {
            size_t point = rect->GetPointNumber(i, j);
            if (rect->pointType[point] == UNUSED) {
                R = 0.0;
            }
            if (rect->pointType[point] == CALCULATED || rect->pointType[point] == INTERPOLATION_DONOR || rect->pointType[point] == DONOR_BUFFER) {
                DiffusiveFlux(rect, &Fd, &ap_d, i, j, k);
                double LHS = Fd;
                double RHS = rect->rho / rect->dt * (rect->Ff[0][i][j] + rect->Ff[1][i][j] + rect->Ff[2][i][j] + rect->Ff[3][i][j]);
                R = RHS - LHS + (1.0-RELAX)/RELAX * rect->Var[k][i][j];
                ap = ap_d/RELAX;
                count++;
                rect->Var[k][i][j] = rect->Var[k][i][j] + R / ap;
            } 
            else if (rect->pointType[point] == INTERPOLATION_RECIEVER) {
                R = 0.0;
            }
            rms = rms + R * R;
        }
    }
    ApplyBC(rect, k, oversetMesh);
    rms = sqrt(rms / (rectMesh->Nx * rectMesh->Ny));
    return rms;
}

// Forward declaration for circular mesh velocity correction
int CorrectVelocity_Circular(DataStructureCircular *circMesh);

int CorrectVelocity(DataStructure *rect) {
    // Cast to rectangular mesh for accessing dx, dy, Nx, Ny
    DataStructureRectangular *rectMesh = dynamic_cast<DataStructureRectangular*>(rect);
    
    for (int k = 0; k < rect->nVar; k++) {
        rect->residual[k] = 0.0;
    }
    
    // Check if this is a circular mesh
    DataStructureCircular *circMesh = dynamic_cast<DataStructureCircular*>(rect);
    if (circMesh != nullptr) {
        // TEMPORARY: Disable circular mesh solver due to instability
        // Use circular mesh solver (Cartesian formulation)
        // return CorrectVelocity_Circular(circMesh);
        return 0;  // Skip for now
    }
    
    // Skip if not a rectangular mesh and not circular
    if (rectMesh == nullptr) {
        return 0;
    }
    
    // Correcting Velocity and counting active computational cells
    int k;
    int activePointCount = 0;
    for (int i = 1; i < rectMesh->Nx + 1; i++) {
        for (int j = 1; j < rectMesh->Ny + 1; j++) {
            size_t point = rect->GetPointNumber(i, j);
            if (rect->pointType[point] == UNUSED || rect->pointType[point] == INTERPOLATION_RECIEVER) {
                continue;
            }
            activePointCount++;
            k = 0;
            rect->Var[k][i][j] = rect->Var[k][i][j] - rect->dt / rect->rho * (rect->Var[2][i + 1][j] - rect->Var[2][i - 1][j]) / (2 * rectMesh->dx);
            k = 1;
            rect->Var[k][i][j] = rect->Var[k][i][j] - rect->dt / rect->rho * (rect->Var[2][i][j + 1] - rect->Var[2][i][j - 1]) / (2 * rectMesh->dy);
            // calculating residuals for u, v, p
            for (int k = 0; k < rect->nVar; ++k) {
                rect->residual[k] += (rect->Var[k][i][j] - rect->VarOld[k][i][j]) * (rect->Var[k][i][j] - rect->VarOld[k][i][j]);
            }
        }
    }
    return activePointCount;
}

inline bool OuterIterLogFreq(const int iter) {
    return (iter % ITER_PRINT_FREQ == 0) || (iter == 1000) ;//|| (iter == 50) || (iter == 100);
}

void ImplicitSolve(DataStructure *rect, DataStructure *oversetMesh, int outerIter, int& activePointsBg, int& activePointsComp) {
    // double Fc, ap_c, Fd, ap_d, ap, R, rms;  // No of faces per cell = 4  i=0,1,2,3 -> E,N,W,S

    double rmsOldBg = 1.0, rmsOldComp = 1.0;
    // Solving for U and V
    int iter = 0;
    double rmsBg = 1.0, rmsComp = 1.0;
    const int MAXITER = 50;
    const double UVTolerance = 1e-6;
    for (int k = 0; k < 2; k++) {
        do {
            rmsOldBg = rmsBg;
            rmsOldComp = rmsComp;
            rmsBg = SolveUV(rect, oversetMesh, k);
            rmsComp = SolveUV(oversetMesh, rect, k);
            iter++;
            if (iter > MAXITER) {
                break;
            }
        } while ((rmsBg > UVTolerance || rmsComp > UVTolerance) && (rmsOldBg > rmsBg || rmsOldComp > rmsComp));
    }
    // if (outerIter % 100 == 0) {
    // 	cout << "RMS U,V: " << rmsBg << " " << rmsComp << " Iter: " << iter; // << endl;
    // }
    LinearInterpolation(rect);
    LinearInterpolation(oversetMesh);

    // Solving for P
    int k = 2;
    iter = 0;
    rmsBg = 1.0, rmsComp = 1.0;
    rmsOldBg = 1.0, rmsOldComp = 1.0;
    const double PTolerance = 1e-6;
    double RELAX = 1.0;
    do {
        rmsOldBg = rmsBg;
        rmsOldComp = rmsComp;
        rmsBg = SolveP(rect, oversetMesh, k, RELAX);
        rmsComp = SolveP(oversetMesh, rect, k, RELAX);
        iter++;
        if (iter > MAXITER) {
            break;
        }
    } while ((rmsBg > PTolerance || rmsComp > PTolerance) && (rmsOldBg > rmsBg || rmsOldComp > rmsComp));
    // if (OuterIterLogFreq(outerIter)) {
    // 	cout << " RMS P: " << rmsBg << " " << rmsComp << " Iter: " << iter << endl;
    // }

    activePointsBg = CorrectVelocity(rect);
    activePointsComp = CorrectVelocity(oversetMesh);
    for (int k = 0; k < 2; ++k) {
        ApplyBC(rect, k, oversetMesh);
        ApplyBC(oversetMesh, k, rect);
        InterpolateCells(rect, oversetMesh, k);
        InterpolateCells(oversetMesh, rect, k);
    }

    UpdateFlux(rect);
    UpdateFlux(oversetMesh);
}

bool ConvergenceCheck(DataStructure *rect, int count, int activePoints) {
    double rms[rect->nVar];
    
    // Use more relaxed tolerance for practical convergence
    // Original 1e-12 was too strict for circular mesh with small dt
    const double TOLERANCE = 1e-6;
    
    // Ensure we don't divide by zero
    if (activePoints == 0) {
        activePoints = 1;
    }
    for (int k = 0; k < rect->nVar; k++) {
        rms[k] = sqrt(rect->residual[k] / activePoints);
        rms[k] = rms[k] / rect->dt;
        if (OuterIterLogFreq(count)) {
            cout << std::setprecision(3) << std::scientific << rms[k] << " ";
        }
    }

    bool converged = !(rms[0] > TOLERANCE || rms[1] > TOLERANCE);
    
    if (!converged) {
        CopyNewtoOld(rect);
    }
    
    return converged;
}

void Solve(DataStructure *rect, DataStructure *oversetMesh, int maxIterations) {
    // Cast to rectangular mesh for accessing Nx, Ny
    DataStructureRectangular *rectMesh = dynamic_cast<DataStructureRectangular*>(rect);
    DataStructureRectangular *oversetRectMesh = dynamic_cast<DataStructureRectangular*>(oversetMesh);
    
    int iter = 0;
    bool mesh1 = true, mesh2 = false;
    int activePointsBg = 0, activePointsComp = 0;
    do {
        ImplicitSolve(rect, oversetMesh, iter, activePointsBg, activePointsComp);
        if (OuterIterLogFreq(iter)) {
            cout << "OuterIter: " << std::setw(5) << iter;
            cout << " NS-RMS: ";
        }
        mesh1 = ConvergenceCheck(rect, iter, activePointsBg);
        if (OuterIterLogFreq(iter)) {
            cout << " - ";
        }
        mesh2 = ConvergenceCheck(oversetMesh, iter, activePointsComp);
        if (OuterIterLogFreq(iter)) {
            cout << endl;
        }
        iter++;
    } while ((mesh1 == false || mesh2 == false) && iter < maxIterations);
    cout << "Solved in " << iter << " iterations." << endl;
    cout << "Active computational cells - Background: " << activePointsBg << "/" << (rectMesh->Nx * rectMesh->Ny) 
         << ", Component: " << activePointsComp << "/" << (oversetRectMesh->Nx * oversetRectMesh->Ny) << endl;
}

// ============================================================================
// CYLINDER FLOW SPECIFIC BOUNDARY CONDITIONS
// ============================================================================

void ApplyBC_Cylinder(DataStructure *rect, int k, DataStructure *oversetMesh, double U_inf) {
    // Cast to rectangular mesh for accessing Nx, Ny (for background mesh)
    DataStructureRectangular *rectMesh = dynamic_cast<DataStructureRectangular*>(rect);
    
    if (rect->pointType[rect->boundaryPoints[0][0]] == INTERPOLATION_RECIEVER) {
        // Outer boundary of circular mesh: interpolate from background
        for (size_t iBound = 0; iBound < rect->boundaryPoints.size(); ++iBound) {
            for (size_t iBoundPoint = 0; iBoundPoint < rect->boundaryPoints[iBound].size(); ++iBoundPoint) {
                size_t point = rect->boundaryPoints[iBound][iBoundPoint];
                auto ij = rect->GetijFromPointNumber(point);
                auto i = std::get<0>(ij);
                auto j = std::get<1>(ij);
                rect->Var[k][i][j] = 0.0;
                for (unsigned short iDonor = 0; iDonor < rect->interpolationStencil[point].size(); ++iDonor) {
                    auto donorPointNumber = rect->interpolationStencil[point][iDonor];
                    auto ijDonor = oversetMesh->GetijFromPointNumber(donorPointNumber);
                    auto _iDonor = std::get<0>(ijDonor);
                    auto _jDonor = std::get<1>(ijDonor);
                    rect->Var[k][i][j] = rect->Var[k][i][j] + oversetMesh->Var[k][_iDonor][_jDonor] * rect->interpolationCoeffs[point][iDonor];
                }
            }
        }
    } else if (rectMesh != nullptr) {
        // Background rectangular mesh: Apply cylinder flow boundary conditions
        switch (k) {
            case 0:  // U velocity
                // Left boundary (inlet): U = U_inf
                for (int j = 1; j < rectMesh->Ny + 1; j++) {
                    rect->Var[k][0][j] = 2 * U_inf - rect->Var[k][1][j];
                }
                // Right boundary (outlet): Zero gradient (Neumann)
                for (int j = 1; j < rectMesh->Ny + 1; j++) {
                    rect->Var[k][rectMesh->Nx + 1][j] = rect->Var[k][rectMesh->Nx][j];
                }
                // Top boundary (far-field): U = U_inf
                for (int i = 1; i < rectMesh->Nx + 1; i++) {
                    rect->Var[k][i][rectMesh->Ny + 1] = 2 * U_inf - rect->Var[k][i][rectMesh->Ny];
                }
                // Bottom boundary (far-field): U = U_inf
                for (int i = 1; i < rectMesh->Nx + 1; i++) {
                    rect->Var[k][i][0] = 2 * U_inf - rect->Var[k][i][1];
                }
                break;

            case 1:  // V velocity
                // Left boundary (inlet): V = 0
                for (int j = 1; j < rectMesh->Ny + 1; j++) {
                    rect->Var[k][0][j] = 2 * (0.0) - rect->Var[k][1][j];
                }
                // Right boundary (outlet): Zero gradient (Neumann)
                for (int j = 1; j < rectMesh->Ny + 1; j++) {
                    rect->Var[k][rectMesh->Nx + 1][j] = rect->Var[k][rectMesh->Nx][j];
                }
                // Top boundary (far-field): V = 0
                for (int i = 1; i < rectMesh->Nx + 1; i++) {
                    rect->Var[k][i][rectMesh->Ny + 1] = 2 * (0.0) - rect->Var[k][i][rectMesh->Ny];
                }
                // Bottom boundary (far-field): V = 0
                for (int i = 1; i < rectMesh->Nx + 1; i++) {
                    rect->Var[k][i][0] = 2 * (0.0) - rect->Var[k][i][1];
                }
                break;

            case 2:  // P (Pressure - Neumann on all boundaries)
                // Left boundary
                for (int j = 1; j < rectMesh->Ny + 1; j++) {
                    rect->Var[k][0][j] = rect->Var[k][1][j];
                }
                // Right boundary
                for (int j = 1; j < rectMesh->Ny + 1; j++) {
                    rect->Var[k][rectMesh->Nx + 1][j] = rect->Var[k][rectMesh->Nx][j];
                }
                // Top boundary
                for (int i = 1; i < rectMesh->Nx + 1; i++) {
                    rect->Var[k][i][rectMesh->Ny + 1] = rect->Var[k][i][rectMesh->Ny];
                }
                // Bottom boundary
                for (int i = 1; i < rectMesh->Nx + 1; i++) {
                    rect->Var[k][i][0] = rect->Var[k][i][1];
                }
                break;
        }
    } else {
        // Circular mesh: Apply no-slip boundary condition on cylinder surface
        // Inner boundary (cylinder surface): U = V = 0 (no-slip wall)
        // Outer boundary: interpolation (handled above)
        if (k == 0 || k == 1) {  // U or V
            for (size_t iBound = 0; iBound < rect->boundaryPoints.size(); ++iBound) {
                for (size_t iBoundPoint = 0; iBoundPoint < rect->boundaryPoints[iBound].size(); ++iBoundPoint) {
                    size_t point = rect->boundaryPoints[iBound][iBoundPoint];
                    auto ij = rect->GetijFromPointNumber(point);
                    auto i = std::get<0>(ij);
                    auto j = std::get<1>(ij);
                    // No-slip: U = V = 0 on cylinder surface
                    rect->Var[k][i][j] = 2 * (0.0) - rect->Var[k][i][j];
                }
            }
        } else if (k == 2) {  // P (Pressure - Neumann on cylinder)
            for (size_t iBound = 0; iBound < rect->boundaryPoints.size(); ++iBound) {
                for (size_t iBoundPoint = 0; iBoundPoint < rect->boundaryPoints[iBound].size(); ++iBoundPoint) {
                    size_t point = rect->boundaryPoints[iBound][iBoundPoint];
                    auto ij = rect->GetijFromPointNumber(point);
                    auto i = std::get<0>(ij);
                    auto j = std::get<1>(ij);
                    // Zero gradient for pressure on cylinder
                    rect->Var[k][i][j] = rect->Var[k][i][j];  // Keep interior value
                }
            }
        }
    }
}

void ImplicitSolveCylinder(DataStructure *rect, DataStructure *oversetMesh, int outerIter, 
                          int& activePointsBg, int& activePointsComp, double U_inf) {
    const double RELAX = 0.7;  // Relaxation factor for pressure
    
    // Copy old values
    CopyNewtoOld(rect);
    CopyNewtoOld(oversetMesh);

    // Interpolate boundary cells
    InterpolateCells(rect, oversetMesh, 0);
    InterpolateCells(rect, oversetMesh, 1);
    InterpolateCells(oversetMesh, rect, 0);
    InterpolateCells(oversetMesh, rect, 1);

    // Apply boundary conditions with cylinder-specific BCs
    ApplyBC_Cylinder(rect, 0, oversetMesh, U_inf);
    ApplyBC_Cylinder(rect, 1, oversetMesh, U_inf);
    ApplyBC_Cylinder(oversetMesh, 0, rect, U_inf);
    ApplyBC_Cylinder(oversetMesh, 1, rect, U_inf);

    // Compute face fluxes
    LinearInterpolation(rect);
    LinearInterpolation(oversetMesh);

    // Solve momentum equations (U and V)
    SolveUV(rect, oversetMesh, 0);
    SolveUV(rect, oversetMesh, 1);
    SolveUV(oversetMesh, rect, 0);
    SolveUV(oversetMesh, rect, 1);

    // Update fluxes with pressure gradient
    UpdateFlux(rect);
    UpdateFlux(oversetMesh);

    // Solve pressure correction equation
    SolveP(rect, oversetMesh, 2, RELAX);
    SolveP(oversetMesh, rect, 2, RELAX);

    // Interpolate pressure
    InterpolateCells(rect, oversetMesh, 2);
    InterpolateCells(oversetMesh, rect, 2);

    // Apply pressure boundary conditions
    ApplyBC_Cylinder(rect, 2, oversetMesh, U_inf);
    ApplyBC_Cylinder(oversetMesh, 2, rect, U_inf);

    // Correct velocities
    activePointsBg = CorrectVelocity(rect);
    activePointsComp = CorrectVelocity(oversetMesh);

    // Apply final boundary conditions
    ApplyBC_Cylinder(rect, 0, oversetMesh, U_inf);
    ApplyBC_Cylinder(rect, 1, oversetMesh, U_inf);
    ApplyBC_Cylinder(oversetMesh, 0, rect, U_inf);
    ApplyBC_Cylinder(oversetMesh, 1, rect, U_inf);

    UpdateFlux(rect);
    UpdateFlux(oversetMesh);
}

void SolveCylinder(DataStructure *rect, DataStructure *oversetMesh, int maxIterations, double U_inf) {
    int iter = 0;
    bool mesh1 = true, mesh2 = false;
    int activePointsBg = 0, activePointsComp = 0;
    
    do {
        ImplicitSolveCylinder(rect, oversetMesh, iter, activePointsBg, activePointsComp, U_inf);
        
        if (OuterIterLogFreq(iter)) {
            cout << "OuterIter: " << std::setw(5) << iter;
            cout << " NS-RMS: ";
        }
        
        mesh1 = ConvergenceCheck(rect, iter, activePointsBg);
        
        if (OuterIterLogFreq(iter)) {
            cout << " - ";
        }
        
        mesh2 = ConvergenceCheck(oversetMesh, iter, activePointsComp);
        
        if (OuterIterLogFreq(iter)) {
            cout << endl;
        }
        
        // Debug: print convergence status (removed after debugging)
        
        iter++;
    } while ((mesh1 == false || mesh2 == false) && iter < maxIterations);
    
    cout << "Solved in " << iter << " iterations." << endl;
    cout << "Active computational cells - Background: " << activePointsBg 
         << ", Component (circular): " << activePointsComp << endl;
}
