//*****************************************************************************
// Title		:src/cpp/main.cpp
// Author	:Tanabe Yuta
// Date		:2019/10/25
// Copyright	:(C)2019 TanabeYuta
//*****************************************************************************

#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>

#include "../src/EllipticalMesher.h"
#include "../src/Point.h"

using namespace EllipticalPAN;
int main() {
    const double height = 2.0, width = 10.0, lead_offset = 0.5;
    const int trail_num = 70, thick_num = 30, lead_num = 77;

    std::ifstream ifs("sample/NACA4412.dat");
    if (!ifs) {
        std::exit(1);
    }
    std::vector<Point<double>> FoilPoints(0);
    while (!ifs.eof()) {
        double x, y;
        ifs >> x >> y;
        FoilPoints.push_back(Point<double>(x, y));
    }

    int foil_num = FoilPoints.size(),
        side_num = trail_num + (foil_num - lead_num) / 2;

    EllipticalMesher<double> Mesh(foil_num + trail_num * 2, thick_num, 1.0e-8);

    Point<double> P0(width, 0.0);
    Point<double> P1(width, height);
    Point<double> P2(lead_offset, height);
    Point<double> P3(lead_offset, -height);
    Point<double> P4(width, -height);

    int ipb = 0;
    for (int i = 0; i < trail_num; i++, ipb++) {
        Mesh.SetBoundaryPoint(ipb, P0 + (FoilPoints[foil_num - 1] - P0) *
                                            (double)(i / (double)trail_num));
    }
    for (int i = foil_num - 1; i >= 0; i--, ipb++) {
        Mesh.SetBoundaryPoint(ipb, FoilPoints[i]);
    }
    for (int i = 0; i < trail_num - 1; i++, ipb++) {
        Mesh.SetBoundaryPoint(
            ipb, FoilPoints[0] + (P0 - FoilPoints[0]) *
                                     (double)((i + 1) / (double)trail_num));
    }
    for (int i = 0; i < thick_num - 1; i++, ipb++) {
        Mesh.SetBoundaryPoint(
            ipb, P0 + (P1 - P0) * (double)(i / (double)(thick_num - 1)));
    }
    for (int i = 0; i < side_num; i++, ipb++) {
        Mesh.SetBoundaryPoint(ipb,
                              P1 + (P2 - P1) * (double)(i / (double)side_num));
    }
    for (int i = 0; i < lead_num; i++, ipb++) {
        Mesh.SetBoundaryPoint(
            ipb,
            Point<double>(
                lead_offset -
                    height * sin(M_PI * (double)(i / (double)(lead_num - 1))),
                height * cos(M_PI * (double)(i / (double)(lead_num - 1)))));
    }
    for (int i = 0; i < side_num - 1; i++, ipb++) {
        Mesh.SetBoundaryPoint(
            ipb, P3 + (P4 - P3) * (double)((i + 1) / (double)side_num));
    }
    for (int i = 0; i < thick_num - 1; i++, ipb++) {
        Mesh.SetBoundaryPoint(
            ipb, P4 + (P0 - P4) * (double)(i / (double)(thick_num - 1)));
    }
    Mesh.MakeMesh();
    Mesh.ExportToVTK("sample/mesh");
    return 0;
}
