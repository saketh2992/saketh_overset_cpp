#include "io/output.h"
#include <cstdio>
#include <fstream>
#include <string>

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

static void write_xyz_block(const DataStructure &m, int k, const std::string &fname) {
    // Writes x y value for all active points (pointType != UNUSED)
    std::ofstream ofs(fname);
    if (!ofs.is_open()) return;
    const int Nx = m.Nx + 2;
    const int Ny = m.Ny + 2;
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            const size_t p = m.GetPointNumber(i, j);
            if (m.pointType[p] == UNUSED) continue;
            double x = m.pointCoordinates[0][p];
            double y = m.pointCoordinates[1][p];
            double v = m.Var[k][i][j];
            ofs << x << " " << y << " " << v << "\n";
        }
        ofs << "\n"; // blank line to separate rows
    }
}

void WriteGnuplotAll(const DataStructure &bg,
                     const DataStructure &comp,
                     const std::string &title,
                     const std::string &outfilePrefix,
                     const std::string &outputDir) {
    // 1) Write XYZ data files for each variable and mesh
    const char* varNames[3] = {"U", "V", "P"};
    for (int k = 0; k < 3; ++k) {
        write_xyz_block(bg, k, outputDir + "/bg_" + varNames[k] + ".xyz");
        write_xyz_block(comp, k, outputDir + "/comp_" + varNames[k] + ".xyz");
    }

    // 2) Write a gnuplot script that builds contours into table files then overlays them in 2D
    std::ofstream plt(outputDir + "/plot_contours.plt");
    if (!plt.is_open()) return;
    plt << "# Auto-generated gnuplot script\n";
    plt << "set term pngcairo size 1200,1000 noenhanced\n";
    plt << "set grid\n";
    plt << "set key outside right\n";
    plt << "set view map\n";
    plt << "unset surface\n";
    plt << "set contour base\n";
    plt << "set cntrparam levels 20\n";
    plt << "set tics out\n";
    plt << "set size ratio -1\n";
    plt << "set style line 1 lc rgb 'black' lw 1.5\n";
    plt << "set style line 2 lc rgb 'red' lw 1.5\n";
    plt << "set dgrid3d 100,100 splines\n";

    for (int k = 0; k < 3; ++k) {
        // Extract contour lines for each mesh separately to table files
        plt << "set table '" << outputDir << "/bg_" << varNames[k] << "_cont.dat'\n";
        plt << "splot '" << outputDir << "/bg_" << varNames[k] << ".xyz' u 1:2:3\n";
        plt << "unset table\n";
        plt << "set table '" << outputDir << "/comp_" << varNames[k] << "_cont.dat'\n";
        plt << "splot '" << outputDir << "/comp_" << varNames[k] << ".xyz' u 1:2:3\n";
        plt << "unset table\n";

        // Plot overlay of contour lines with different colors
        plt << "set output '" << outputDir << "/" << outfilePrefix << "_" << varNames[k] << ".png'\n";
        plt << "set title '" << title << " - " << varNames[k] << "'\n";
        plt << "plot \\\n+  '" << outputDir << "/bg_" << varNames[k] << "_cont.dat' w l ls 1 title 'bg', \\\n+  '" << outputDir << "/comp_" << varNames[k] << "_cont.dat' w l ls 2 title 'comp'\n";
        plt << "unset output\n";
    }

    plt << "unset dgrid3d\n";
    plt.close();
}
