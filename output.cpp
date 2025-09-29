#include "output.h"
#include <cstdio>

void GetOutput(DataStructure *rect, std::string input) {
    FILE *fp;
    fp = fopen(input.c_str(), "w");

    for (int k = 0; k < rect->nVar; k++) {
        fprintf(fp, "\n ########## Data for k = %d ############ \n", k);
        for (int j = 0; j < rect->Ny + 2; j++) {
            for (int i = 0; i < rect->Nx + 2; i++) {
                fprintf(fp, "%lf \t", rect->Var[k][i][j]);
            }
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
}