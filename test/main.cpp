//*****************************************************************************
//Title		:src/cpp/main.cpp
//Author	:Tanabe Yuta
//Date		:2019/10/25
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "EllipticalMesher.h"
#include "Point.h"


using namespace EllipticalPAN;
int main() {
	EllipticalMesher<double> Mesh(50, 199, 1.0e-8);		

	Point<double> P0(0.0, 0.0);
	Point<double> P1(4.0, 0.0);
	Point<double> P2(4.0, 4.0);
	Point<double> P3(6.0, 6.0);
	Point<double> P4(10.0, 6.0);
	Point<double> P5(10.0, 10.0);
	Point<double> P6(6.0, 10.0);
	Point<double> P7(0.0, 4.0);

	int ipb = 0;
	for (int i = 0; i < 49; i++, ipb++) {
		Mesh.SetBoundaryPoint(ipb, P0 + (P1 - P0)*i / 49);
	}
	for (int i = 0; i < 49; i++, ipb++) {
		Mesh.SetBoundaryPoint(ipb, P1 + (P2 - P1)*i / 49);
	}
	for (int i = 0; i < 100; i++, ipb++) {
		Mesh.SetBoundaryPoint(ipb, P2 + Point<double>(2.0*(1.0 - cos(0.5*M_PI*i / 100)), 2.0*sin(0.5*M_PI*i / 100)));
	}
	for (int i = 0; i < 49; i++, ipb++) {
		Mesh.SetBoundaryPoint(ipb, P3 + (P4 - P3)*i / 49);
	}
	for (int i = 0; i < 49; i++, ipb++) {
		Mesh.SetBoundaryPoint(ipb, P4 + (P5 - P4)*i / 49);
	}
	for (int i = 0; i < 49; i++, ipb++) {
		Mesh.SetBoundaryPoint(ipb, P5 + (P6 - P5)*i / 49);
	}
	for (int i = 0; i < 100; i++, ipb++) {
		Mesh.SetBoundaryPoint(ipb, P6 + Point<double>(-6.0*sin(0.5*M_PI*i / 100), -6.0*(1.0 - cos(0.5*M_PI*i / 100))));
	}
	for (int i = 0; i < 49; i++, ipb++) {
		Mesh.SetBoundaryPoint(ipb, P7 + (P0 - P7)*i / 49);
	}

	Mesh.MakeMesh();					
	Mesh.ExportToVTK("sample/mesh");			
	Mesh.ExportForPANSFEM("sample/mesh");
	return 0;
}