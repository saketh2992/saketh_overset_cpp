#ifndef OUTPUT_H
#define OUTPUT_H

#include <string>
#include "mesh/datastructure.h"

// Output functions
void GetOutput(DataStructure *rect, std::string input);

// Post-processing: write gnuplot inputs and a plotting script for contours
void WriteGnuplotAll(const DataStructure &bg,
					 const DataStructure &comp,
					 const std::string &title = "NS contour for 2d plate",
					 const std::string &outfilePrefix = "lid_overset",
					 const std::string &outputDir = "TEMP");

#endif // OUTPUT_H