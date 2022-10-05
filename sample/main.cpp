#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "../src/mesher.h"
#include "../src/vec2.h"

using namespace EllipticalPAN;

int main() {
    const double height = 2.0, width = 10.0, offset = 0.5;
    const int trail_num = 70, thick_num = 30, lead_num = 77;

    std::ifstream ifs("NACA4412.dat");
    if (!ifs) {
        std::exit(1);
    }
    std::vector<Vec2<double>> PF(0);
    while (!ifs.eof()) {
        double x, y;
        ifs >> x >> y;
        PF.push_back(Vec2<double>(x, y));
    }

    int foil_num = PF.size(), side_num = trail_num + (foil_num - lead_num) / 2;

    Mesher<Vec2<double>, double> Mesh(foil_num + trail_num * 2, thick_num);

    Vec2<double> P0(width, 0.0);
    Vec2<double> P1(width, height);
    Vec2<double> P2(offset, height);
    Vec2<double> P3(offset, -height);
    Vec2<double> P4(width, -height);

    int ipb = 0;
    for (int i = 0; i < trail_num; i++, ipb++) {
        Mesh.SetPoint(
            ipb, 0,
            P0 + (PF[foil_num - 1] - P0) * (double)(i / (double)trail_num));
    }
    for (int i = foil_num - 1; i >= 0; i--, ipb++) {
        Mesh.SetPoint(ipb, 0, PF[i]);
    }
    for (int i = 0; i < trail_num; i++, ipb++) {
        Mesh.SetPoint(
            ipb, 0,
            PF[0] + (P0 - PF[0]) * (double)((i + 1) / (double)trail_num));
    }
    ipb = 0;
    for (int i = 0; i < thick_num; i++, ipb++) {
        Mesh.SetPoint(-1, ipb,
                      P0 + (P1 - P0) * (double)(i / (double)(thick_num - 1)));
    }
    ipb = 0;
    for (int i = 0; i < side_num; i++, ipb++) {
        Mesh.SetPoint(ipb, -1, P4 + (P3 - P4) * (double)(i / (double)side_num));
    }
    for (int i = 0; i < lead_num; i++, ipb++) {
        Mesh.SetPoint(
            ipb, -1,
            Vec2<double>(
                offset - height * sin(M_PI * i / (double)(lead_num - 1)),
                -height * cos(M_PI * i / (double)(lead_num - 1))));
    }
    for (int i = 0; i < side_num; i++, ipb++) {
        Mesh.SetPoint(ipb, -1,
                      P2 + (P1 - P2) * (double)((i + 1) / (double)side_num));
    }
    ipb = 0;
    for (int i = 0; i < thick_num; i++, ipb++) {
        Mesh.SetPoint(0, ipb,
                      P0 + (P4 - P0) * (double)(i / (double)(thick_num - 1)));
    }

    std::cout << "Converged at " << Mesh.Generate() << std::endl;

    std::ofstream fout("mesh.vtk");
    Mesh.ExportAsVTK(fout);
    fout.close();
    return 0;
}
