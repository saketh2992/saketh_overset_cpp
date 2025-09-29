#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <limits>
#include <cstddef>

// Type definitions
using su2double = double;
using std::size_t; // ensure size_t is visible without std:: prefix

// Constants
const size_t ITER_PRINT_FREQ = 5000;
const double PI = 3.141592653589793238463;
const size_t NS_MAXITER = 20;
const double TOLERANCE = 1024 * std::numeric_limits<double>::epsilon();

// SU2 specific constants
const size_t SU2_BBOX_SIZE = 4; /*!< \brief Size of the bounding box array that is allocated for each element*/
constexpr su2double su2double_lowest = std::numeric_limits<su2double>::lowest();
constexpr su2double su2double_highest = std::numeric_limits<su2double>::max();

// Point type enumeration
enum POINT_TYPE {
    UNUSED = 0,
    CALCULATED = 1,
    INTERPOLATION_DONOR = 2,
    INTERPOLATION_RECIEVER = 3,
    DONOR_BUFFER = 4,
    BC_SPECIFIED = 5,
};

#endif // CONSTANTS_H