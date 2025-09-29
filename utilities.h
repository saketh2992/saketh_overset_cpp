#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include "constants.h"

// Mathematical utility functions
std::vector<double> solveAxB(const std::vector<std::vector<double>>& A, const std::vector<double>& B);
std::vector<double> matAdd(const std::vector<double>& A, const std::vector<double>& B);
std::vector<double> interpolation_coeff(const std::vector<std::vector<double>>& points, const std::vector<double>& xy);

#endif // UTILITIES_H