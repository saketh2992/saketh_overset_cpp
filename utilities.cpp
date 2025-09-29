#include "utilities.h"
#include <cmath>
#include <limits>

using namespace std;

vector<double> solveAxB(const vector<vector<double> > &A, const vector<double> &B) {
    size_t n = A.size();
    vector<double> x(n, 0.0);

    // Pre-calculate inverses of diagonal elements (avoid division within the loop)
    vector<double> inv_diag(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        inv_diag[i] = 1.0 / A[i][i];
    }

    // Loop for iterations, Guass Seidel
    size_t GS_itr = 0;
    for (; GS_itr < NS_MAXITER; ++GS_itr) {
        // Update each element of the solution vector
        double sum = 0.0;
        bool converged = true;
        for (size_t j = 0; j < n; ++j) {
            sum = 0.0;
            for (size_t k = 0; k < n; ++k) {
                if (j != k) {
                    sum += A[j][k] * x[k];
                }
            }
            x[j] = (B[j] - sum) * inv_diag[j];
            converged &= (abs(x[j] * A[j][j] + sum - B[j]) <= TOLERANCE);
        }

        if (converged) {
            break;
        }
    }
    // cout << "GS=" << GS_itr << ", ";
    return x;
}

vector<double> matAdd(const vector<double> &A, const vector<double> &B) {
    size_t m = A.size();
    vector<double> C(m, 0.0);
    for (size_t i = 0; i < m; i++) {
        C[i] = A[i] + B[i];
    }
    return C;
}

vector<double> interpolation_coeff(const vector<vector<double> >& points, const vector<double>& xy) {
    auto transform = [](const vector<double> &coeff, const vector<double> &points) -> double {
        const double phi = coeff[0], chi = coeff[1];
        return (1.0 - phi) * (1.0 - chi) * points[0] + (phi) * (1.0 - chi) * points[1] + (phi) * (chi)*points[2] + (1.0 - phi) * (chi)*points[3];
    };
    auto partialDerivative_phi = [](const vector<double> &coeff, const vector<double> &points) {
        const double chi = coeff[1];
        return (-1.0) * (1.0 - chi) * points[0] + (1.0 - chi) * points[1] + (chi)*points[2] + (-1.0) * (chi)*points[3];
    };
    auto partialDerivative_chi = [](const vector<double> &coeff, const vector<double> &points) {
        const double phi = coeff[0];
        return (1.0 - phi) * (-1.0) * points[0] + (phi) * (-1.0) * points[1] + (phi)*points[2] + (1.0 - phi) * points[3];
    };

    vector<double> coeff = {0.5, 0.5};  // initial guess value of phi,

    // NR, Newton raphson to minimze B, Ax = B
    size_t NR_itr = 0;
    for (; NR_itr < NS_MAXITER; NR_itr++) {
        vector<vector<double> > A = {{partialDerivative_phi(coeff, points[0]), partialDerivative_chi(coeff, points[0])},
                                     {partialDerivative_phi(coeff, points[1]), partialDerivative_chi(coeff, points[1])}};
        vector<double> B = {-1.0 * transform(coeff, points[0]) + xy[0], -1.0 * transform(coeff, points[1]) + xy[1]};
        
        coeff = matAdd(coeff, solveAxB(A, B));
        // check for convergence
        bool converged = true;
        for (int eq = 0; eq < 2; ++eq) {
            converged &= (abs(B[eq]) <= TOLERANCE);
        }
        if (converged) {
            break;
        }
    }
    // cout << "NR=" << NR_itr << ", ";
    return vector<double>{(1.0 - coeff[0]) * (1.0 - coeff[1]), (coeff[0]) * (1.0 - coeff[1]), (coeff[0]) * (coeff[1]), (1.0 - coeff[0]) * (coeff[1])};
}