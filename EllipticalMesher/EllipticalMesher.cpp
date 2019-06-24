//*****************************************************************************
//Title		:EllipticalMesher/EllipticalMesher.cpp
//Author	:Tanabe Yuta
//Date		:2019/04/16
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#include <string>
#include <fstream>
#include <sstream>
#include "EllipticalMesher.h"


EllipticalMesher::EllipticalMesher() : inum(0), jnum(0) {}


EllipticalMesher::~EllipticalMesher(){}


EllipticalMesher::EllipticalMesher(int _inum, int _jnum) : inum(_inum), jnum(_jnum){
	this->pb = std::vector<Point>(2 * (this->inum + this->jnum - 2));
	this->pin = std::vector<std::vector<Point>>(this->inum, std::vector<Point>(this->jnum, Point(0.0, 0.0)));
}


void EllipticalMesher::SetBoundaryPoint(int _i, Point _point){
	this->pb[_i] = _point;
}


void EllipticalMesher::MakeMesh(){
	//----------Dirichlet条件の設定----------
	int pbi = 0;		//境界のカウンタ
	//.....下面.....
	for (int i = 0; i < this->inum - 1; i++, pbi++) {
		this->pin[i][0] = this->pb[pbi];
	}
	//.....右面.....
	for (int i = 0; i < this->jnum - 1; i++, pbi++) {
		this->pin[this->inum - 1][i] = this->pb[pbi];
	}
	//.....上面.....
	for (int i = this->inum - 1; i > 0; i--, pbi++) {
		this->pin[i][this->jnum - 1] = this->pb[pbi];
	}
	//.....左面.....
	for (int i = this->jnum - 1; i > 0; i--, pbi++) {
		this->pin[0][i] = this->pb[pbi];
	}

	//----------Laplace方程式の求解----------
	for (int k = 0; k < this->ITRMAX; k++) {
		double dxmax = 0.0;
		std::cout << "iteration = " << k + 1;

		for (int i = 1; i < this->inum - 1; i++) {
			for (int j = 1; j < this->jnum - 1; j++) {
				//.....Laplace方程式の係数を計算.....
				Point xxi = 0.5*(this->pin[i + 1][j] - this->pin[i - 1][j]);
				Point xiota = 0.5*(this->pin[i][j + 1] - this->pin[i][j - 1]);
				double alpha = pow(xiota.x[0], 2.0) + pow(xiota.x[1], 2.0);
				double beta = xxi.x[0] * xiota.x[0] + xxi.x[1] * xiota.x[1];
				double ganma = pow(xxi.x[0], 2.0) + pow(xxi.x[1], 2.0);

				//.....座標を更新.....
				Point tmp = this->pin[i][j];
				this->pin[i][j] = 0.5*(alpha*(this->pin[i + 1][j] + this->pin[i - 1][j]) - 0.5*beta*(this->pin[i + 1][j + 1] - this->pin[i - 1][j + 1] - this->pin[i + 1][j - 1] + this->pin[i - 1][j - 1]) + ganma * (this->pin[i][j + 1] + this->pin[i][j - 1])) / (alpha + ganma);
				
				//.....座標増分のノルムの最大値を更新.....
				double dxtmp = (tmp - this->pin[i][j]).Norm();
				if (dxmax < dxtmp) {
					dxmax = dxtmp;
				}
			}
		}
		std::cout << "\t" << dxmax << "\n";

		//.....収束判定.....
		if (dxmax < DEPS) {
			std::cout << "Converged:k = " << k + 1;
			break;
		}
	}
}


void EllipticalMesher::ExportToVTK(std::string _fname){
	std::ofstream fout(_fname + ".vtk");
	
	//----------Headerの出力----------
	fout << "# vtk DataFile Version 4.1\n";
	fout << "vtk output\n";
	fout << "ASCII\n";
	fout << "DATASET UNSTRUCTURED_GRID\n";

	//----------点の追加----------
	fout << "\nPOINTS\t" << this->inum * this->jnum << "\tfloat\n";
	for (auto pinrow : this->pin) {
		for (auto pinone : pinrow) {
			fout << pinone.x[0] << "\t" << pinone.x[1] << "\t" << 0.0 << "\n";
		}
	}

	//----------要素の追加----------
	fout << "\nCELLS " << (this->inum - 1)*(this->jnum - 1) << "\t" << (this->inum - 1)*(this->jnum - 1) * 5 << "\n";
	for (int i = 0; i < this->inum - 1; i++) {
		for (int j = 0; j < this->jnum - 1; j++) {
			fout << "4 " << i * (this->jnum)+j << "\t" << i * (this->jnum)+j + 1 << "\t" << (i + 1)*(this->jnum)+j + 1 << "\t" << (i + 1)*(this->jnum)+j << "\n";
		}
	}

	//----------要素タイプの設定----------
	fout << "\nCELL_TYPES\t" << (this->inum - 1)*(this->jnum - 1) << "\n";
	for (int i = 0; i < (this->inum - 1)*(this->jnum - 1); i++) {
		fout << "9\n";
	}

	fout.close();
}


void EllipticalMesher::ExportForPANSFEM(std::string _fnameheader){
	//----------Nodeの出力----------
	std::ofstream fout_node(_fnameheader + "_Node.csv");

	//.....Headerの出力.....
	fout_node << "ID,DOX,x0,x1";

	//.....座標を出力.....
	for (int i = 0; i < this->inum; i++) {
		for (int j = 0; j < this->jnum; j++) {
			fout_node << std::endl << i * (this->jnum) + j << "," << "2" << "," << this->pin[i][j].x[0] << "," << this->pin[i][j].x[1];
		}
	}

	fout_node.close();

	//----------Elementの出力----------
	std::ofstream fout_element(_fnameheader + "_Element.csv");

	//.....Headerの出力.....
	fout_element << "ID,NON,n0,n1,n2,n3";

	//.....節点番号を出力.....
	for (int i = 0; i < this->inum - 1; i++) {
		for (int j = 0; j < this->jnum - 1; j++) {
			fout_element << std::endl << i * (this->jnum - 1) + j << "," << "4" << "," << i * (this->jnum) + j << "," << (i + 1)*(this->jnum) + j << "," << (i + 1)*(this->jnum) + j + 1 << "," << i * (this->jnum) + j + 1;
		}
	}

	fout_element.close();
}