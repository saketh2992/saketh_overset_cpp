// ============================================================================
// CIRCULAR MESH SOLVER - CARTESIAN FORMULATION
// ============================================================================
// These functions solve the Navier-Stokes equations on the circular mesh
// using Cartesian coordinates (x, y) instead of polar coordinates (r, theta).
// The circular mesh points already have Cartesian coordinates stored in
// pointCoordinates[0] (x) and pointCoordinates[1] (y).
// ============================================================================

#include "mesh/datastructure_circular.h"
#include "util/constants.h"
#include <cmath>
#include <iostream>

// Compute local cell spacing for circular mesh using Cartesian coordinates
void ComputeCircularCellSpacing(DataStructureCircular *circMesh, int i, int j, 
                                double &dx_local, double &dy_local) {
    // Get current point
    size_t p = circMesh->GetPointNumber(i, j);
    double x_center = circMesh->pointCoordinates[0][p];
    double y_center = circMesh->pointCoordinates[1][p];
    
    // Compute spacing in x-direction (between i-1 and i+1)
    if (i > 0 && i < circMesh->Nr + 1) {
        size_t p_left = circMesh->GetPointNumber(i - 1, j);
        size_t p_right = circMesh->GetPointNumber(i + 1, j);
        dx_local = (circMesh->pointCoordinates[0][p_right] - circMesh->pointCoordinates[0][p_left]) / 2.0;
    } else {
        dx_local = circMesh->dr;  // Fallback to radial spacing
    }
    
    // Compute spacing in y-direction (between j-1 and j+1)
    if (j > 0 && j < circMesh->Ntheta + 1) {
        size_t p_down = circMesh->GetPointNumber(i, j - 1);
        size_t p_up = circMesh->GetPointNumber(i, j + 1);
        dy_local = (circMesh->pointCoordinates[1][p_up] - circMesh->pointCoordinates[1][p_down]) / 2.0;
    } else {
        // For periodic boundary in theta direction
        size_t p_down = circMesh->GetPointNumber(i, j - 1 + circMesh->Ntheta);
        size_t p_up = circMesh->GetPointNumber(i, (j + 1) % circMesh->Ntheta);
        dy_local = (circMesh->pointCoordinates[1][p_up] - circMesh->pointCoordinates[1][p_down]) / 2.0;
    }
    
    // Ensure positive spacing
    dx_local = fabs(dx_local);
    dy_local = fabs(dy_local);
    if (dx_local < 1e-10) dx_local = circMesh->dr;
    if (dy_local < 1e-10) dy_local = circMesh->dr * circMesh->dtheta;
}

// Solve momentum equations (U, V) for circular mesh
double SolveUV_Circular(DataStructureCircular *circMesh, DataStructure *oversetMesh, int k) {
    /* k = 0 -> U, k = 1 -> V */
    
    int count = 0;
    double rms = 0.0;
    
    // Loop over interior points (excluding ghost cells)
    for (int i = 1; i <= circMesh->Nr; i++) {
        for (int j = 1; j <= circMesh->Ntheta; j++) {
            size_t point = circMesh->GetPointNumber(i, j);
            
            // Skip unused or interpolation receiver points
            if (circMesh->pointType[point] == UNUSED || 
                circMesh->pointType[point] == INTERPOLATION_RECIEVER) {
                continue;
            }
            
            // Only solve for calculated, donor, and donor buffer points
            if (circMesh->pointType[point] == CALCULATED || 
                circMesh->pointType[point] == INTERPOLATION_DONOR || 
                circMesh->pointType[point] == DONOR_BUFFER) {
                
                // Get local cell spacing
                double dx_local, dy_local;
                ComputeCircularCellSpacing(circMesh, i, j, dx_local, dy_local);
                
                // Neighbor indices with periodic boundary in theta
                int i_next = i + 1;
                int i_prev = i - 1;
                int j_next = (j == circMesh->Ntheta) ? 1 : j + 1;
                int j_prev = (j == 1) ? circMesh->Ntheta : j - 1;
                
                // Convective term (upwind scheme)
                double U_center = circMesh->Var[0][i][j];
                double V_center = circMesh->Var[1][i][j];
                
                double Var_east = circMesh->Var[k][i_next][j];
                double Var_west = circMesh->Var[k][i_prev][j];
                double Var_north = circMesh->Var[k][i][j_next];
                double Var_south = circMesh->Var[k][i][j_prev];
                double Var_center = circMesh->Var[k][i][j];
                
                // Upwind convection
                double conv_x = 0.0, conv_y = 0.0;
                if (U_center > 0) {
                    conv_x = U_center * (Var_center - Var_west) / dx_local;
                } else {
                    conv_x = U_center * (Var_east - Var_center) / dx_local;
                }
                
                if (V_center > 0) {
                    conv_y = V_center * (Var_center - Var_south) / dy_local;
                } else {
                    conv_y = V_center * (Var_north - Var_center) / dy_local;
                }
                
                // Diffusive term (central difference)
                double diff_xx = (Var_east - 2.0 * Var_center + Var_west) / (dx_local * dx_local);
                double diff_yy = (Var_north - 2.0 * Var_center + Var_south) / (dy_local * dy_local);
                double diffusion = circMesh->nu * (diff_xx + diff_yy);
                
                // Pressure gradient (will be added in pressure correction)
                // For now, just solve momentum without pressure gradient
                
                // Time derivative
                double dVar_dt = (Var_center - circMesh->VarOld[k][i][j]) / circMesh->dt;
                
                // Residual: dU/dt + U dU/dx + V dU/dy = nu * (d²U/dx² + d²U/dy²)
                double R = -(dVar_dt + conv_x + conv_y - diffusion);
                
                // Update
                double volume = dx_local * dy_local;
                double ap = volume / circMesh->dt + 
                           fabs(U_center) * volume / dx_local + 
                           fabs(V_center) * volume / dy_local +
                           2.0 * circMesh->nu * volume * (1.0/(dx_local*dx_local) + 1.0/(dy_local*dy_local));
                
                if (ap > 1e-12) {
                    circMesh->Var[k][i][j] = Var_center + 0.7 * R / ap;  // Under-relaxation
                }
                
                count++;
                rms += R * R;
            }
        }
    }
    
    if (count > 0) {
        rms = sqrt(rms / count);
    }
    
    return rms;
}

// Solve pressure Poisson equation for circular mesh
double SolveP_Circular(DataStructureCircular *circMesh, DataStructure *oversetMesh, 
                       int k, double RELAX) {
    /* k = 2 for pressure */
    
    int count = 0;
    double rms = 0.0;
    
    for (int i = 1; i <= circMesh->Nr; i++) {
        for (int j = 1; j <= circMesh->Ntheta; j++) {
            size_t point = circMesh->GetPointNumber(i, j);
            
            if (circMesh->pointType[point] == UNUSED || 
                circMesh->pointType[point] == INTERPOLATION_RECIEVER) {
                continue;
            }
            
            if (circMesh->pointType[point] == CALCULATED || 
                circMesh->pointType[point] == INTERPOLATION_DONOR || 
                circMesh->pointType[point] == DONOR_BUFFER) {
                
                // Get local cell spacing
                double dx_local, dy_local;
                ComputeCircularCellSpacing(circMesh, i, j, dx_local, dy_local);
                
                // Neighbor indices with periodic boundary
                int i_next = i + 1;
                int i_prev = i - 1;
                int j_next = (j == circMesh->Ntheta) ? 1 : j + 1;
                int j_prev = (j == 1) ? circMesh->Ntheta : j - 1;
                
                // Compute velocity divergence
                double U_east = circMesh->Var[0][i_next][j];
                double U_west = circMesh->Var[0][i_prev][j];
                double V_north = circMesh->Var[1][i][j_next];
                double V_south = circMesh->Var[1][i][j_prev];
                
                double div_U = (U_east - U_west) / (2.0 * dx_local) + 
                              (V_north - V_south) / (2.0 * dy_local);
                
                // Source term for pressure correction
                double source = circMesh->rho / circMesh->dt * div_U;
                
                // Pressure Poisson: ∇²P = -ρ/dt * div(U)
                double P_center = circMesh->Var[k][i][j];
                double P_east = circMesh->Var[k][i_next][j];
                double P_west = circMesh->Var[k][i_prev][j];
                double P_north = circMesh->Var[k][i][j_next];
                double P_south = circMesh->Var[k][i][j_prev];
                
                // Laplacian of pressure
                double laplacian_P = (P_east - 2.0 * P_center + P_west) / (dx_local * dx_local) +
                                    (P_north - 2.0 * P_center + P_south) / (dy_local * dy_local);
                
                // Residual
                double R = source - laplacian_P;
                
                // Update pressure with relaxation
                double ap = 2.0 / (dx_local * dx_local) + 2.0 / (dy_local * dy_local);
                if (ap > 1e-12) {
                    circMesh->Var[k][i][j] = P_center + RELAX * R / ap;
                }
                
                count++;
                rms += R * R;
            }
        }
    }
    
    if (count > 0) {
        rms = sqrt(rms / count);
    }
    
    return rms;
}

// Correct velocities using pressure gradient for circular mesh
int CorrectVelocity_Circular(DataStructureCircular *circMesh) {
    
    // Initialize residuals
    for (int k = 0; k < circMesh->nVar; k++) {
        circMesh->residual[k] = 0.0;
    }
    
    int activePointCount = 0;
    
    for (int i = 1; i <= circMesh->Nr; i++) {
        for (int j = 1; j <= circMesh->Ntheta; j++) {
            size_t point = circMesh->GetPointNumber(i, j);
            
            if (circMesh->pointType[point] == UNUSED || 
                circMesh->pointType[point] == INTERPOLATION_RECIEVER) {
                continue;
            }
            
            activePointCount++;
            
            // Get local cell spacing
            double dx_local, dy_local;
            ComputeCircularCellSpacing(circMesh, i, j, dx_local, dy_local);
            
            // Neighbor indices
            int i_next = i + 1;
            int i_prev = i - 1;
            int j_next = (j == circMesh->Ntheta) ? 1 : j + 1;
            int j_prev = (j == 1) ? circMesh->Ntheta : j - 1;
            
            // Pressure gradient
            double dP_dx = (circMesh->Var[2][i_next][j] - circMesh->Var[2][i_prev][j]) / (2.0 * dx_local);
            double dP_dy = (circMesh->Var[2][i][j_next] - circMesh->Var[2][i][j_prev]) / (2.0 * dy_local);
            
            // Correct velocities
            circMesh->Var[0][i][j] -= circMesh->dt / circMesh->rho * dP_dx;
            circMesh->Var[1][i][j] -= circMesh->dt / circMesh->rho * dP_dy;
            
            // Calculate residuals
            for (int k = 0; k < circMesh->nVar; k++) {
                double diff = circMesh->Var[k][i][j] - circMesh->VarOld[k][i][j];
                circMesh->residual[k] += diff * diff;
            }
        }
    }
    
    return activePointCount;
}

// Apply boundary conditions for circular mesh (Cartesian velocities)
void LinearInterpolation_Circular(DataStructureCircular *circMesh) {
    // For circular mesh, we can compute face fluxes using Cartesian coordinates
    // This is a simplified version - just stores velocity at faces
    
    for (int i = 1; i <= circMesh->Nr; i++) {
        for (int j = 1; j <= circMesh->Ntheta; j++) {
            size_t point = circMesh->GetPointNumber(i, j);
            
            if (circMesh->pointType[point] == UNUSED || 
                circMesh->pointType[point] == INTERPOLATION_RECIEVER) {
                continue;
            }
            
            // Simple averaging for face fluxes
            // In a structured circular mesh, we can compute fluxes similarly
            // This is a placeholder - full implementation would compute actual face fluxes
            
            // For now, just ensure the flux arrays are initialized
            if (i < circMesh->Nr && j < circMesh->Ntheta) {
                double U_avg = (circMesh->Var[0][i][j] + circMesh->Var[0][i+1][j]) * 0.5;
                double V_avg = (circMesh->Var[1][i][j] + circMesh->Var[1][i][j+1]) * 0.5;
                
                // Store in flux array if needed
                // circMesh->Ff is defined for face fluxes
            }
        }
    }
}
