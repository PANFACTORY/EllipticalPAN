//*****************************************************************************
//Title		:src/cpp/EllipticalMesher.h
//Author	:Tanabe Yuta
//Date		:2019/10/25
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>


#include "Point.h"


namespace EllipticalPAN{
	template<class T>
	class EllipticalMesher
	{
public:
		EllipticalMesher();
		~EllipticalMesher();
		EllipticalMesher(int _inum, int _jnum, T _DEPS);
		
		
		void SetBoundaryPoint(int, Point<T>);		
		void MakeMesh();						
		void ExportToVTK(std::string);			
		void ExportForPANSFEM(std::string);		


private:
		const int inum;							
		const int jnum;						


		std::vector<Point<T> > pb;					
		std::vector<std::vector<Point<T> > > pin;	


		static const int ITRMAX = 100000;		
		const T DEPS;				
	};


	template<class T>
	EllipticalMesher<T>::EllipticalMesher() : inum(0), jnum(0), DEPS(T()) {}


	template<class T>
	EllipticalMesher<T>::~EllipticalMesher(){}


	template<class T>
	EllipticalMesher<T>::EllipticalMesher(int _inum, int _jnum, T _DEPS) : inum(_inum), jnum(_jnum), DEPS(_DEPS){
		this->pb = std::vector<Point<T> >(2 * (this->inum + this->jnum - 2));
		this->pin = std::vector<std::vector<Point<T> > >(this->inum, std::vector<Point<T> >(this->jnum, Point<T>(T(), T())));
	}


	template<class T>
	void EllipticalMesher<T>::SetBoundaryPoint(int _i, Point<T> _point){
		this->pb[_i] = _point;
	}


	template<class T>
	void EllipticalMesher<T>::MakeMesh(){
		//----------Set Dirichlet condition----------
		int pbi = 0;
		//.....Bottom edge.....
		for (int i = 0; i < this->inum - 1; i++, pbi++) {
			this->pin[i][0] = this->pb[pbi];
		}
		//.....Right edge.....
		for (int i = 0; i < this->jnum - 1; i++, pbi++) {
			this->pin[this->inum - 1][i] = this->pb[pbi];
		}
		//.....Top edge.....
		for (int i = this->inum - 1; i > 0; i--, pbi++) {
			this->pin[i][this->jnum - 1] = this->pb[pbi];
		}
		//.....Left edge.....
		for (int i = this->jnum - 1; i > 0; i--, pbi++) {
			this->pin[0][i] = this->pb[pbi];
		}

		//----------Solve Laplace equation----------
		for (int k = 0; k < this->ITRMAX; k++) {
			T dxmax = T();
			//std::cout << "iteration = " << k + 1;

			for (int i = 1; i < this->inum - 1; i++) {
				for (int j = 1; j < this->jnum - 1; j++) {
					//.....Make Laplace equation.....
					Point<T> xxi = 0.5*(this->pin[i + 1][j] - this->pin[i - 1][j]);
					Point<T> xiota = 0.5*(this->pin[i][j + 1] - this->pin[i][j - 1]);
					T alpha = pow(xiota.x[0], 2.0) + pow(xiota.x[1], 2.0);
					T beta = xxi.x[0] * xiota.x[0] + xxi.x[1] * xiota.x[1];
					T ganma = pow(xxi.x[0], 2.0) + pow(xxi.x[1], 2.0);

					//.....Update values.....
					Point<T> tmp = this->pin[i][j];
					this->pin[i][j] = 0.5*(alpha*(this->pin[i + 1][j] + this->pin[i - 1][j]) - 0.5*beta*(this->pin[i + 1][j + 1] - this->pin[i - 1][j + 1] - this->pin[i + 1][j - 1] + this->pin[i - 1][j - 1]) + ganma * (this->pin[i][j + 1] + this->pin[i][j - 1])) / (alpha + ganma);
					
					//.....Get norm maximam.....
					T dxtmp = (tmp - this->pin[i][j]).Norm();
					if (dxmax < dxtmp) {
						dxmax = dxtmp;
					}
				}
			}
			//std::cout << "\t" << dxmax << "\n";

			//.....Check convergence.....
			if (dxmax < DEPS) {
				std::cout << "Converged:k = " << k + 1;
				break;
			}
		}
	}


	template<class T>
	void EllipticalMesher<T>::ExportToVTK(std::string _fname){
		std::ofstream fout(_fname + ".vtk");
		
		//----------Header�̏o��----------
		fout << "# vtk DataFile Version 4.1\n";
		fout << "vtk output\n";
		fout << "ASCII\n";
		fout << "DATASET UNSTRUCTURED_GRID\n";

		//----------�_�̒ǉ�----------
		fout << "\nPOINTS\t" << this->inum * this->jnum << "\tfloat\n";
		for (auto pinrow : this->pin) {
			for (auto pinone : pinrow) {
				fout << pinone.x[0] << "\t" << pinone.x[1] << "\t" << 0.0 << "\n";
			}
		}

		//----------�v�f�̒ǉ�----------
		fout << "\nCELLS " << (this->inum - 1)*(this->jnum - 1) << "\t" << (this->inum - 1)*(this->jnum - 1) * 5 << "\n";
		for (int i = 0; i < this->inum - 1; i++) {
			for (int j = 0; j < this->jnum - 1; j++) {
				fout << "4 " << i * (this->jnum)+j << "\t" << i * (this->jnum)+j + 1 << "\t" << (i + 1)*(this->jnum)+j + 1 << "\t" << (i + 1)*(this->jnum)+j << "\n";
			}
		}

		//----------�v�f�^�C�v�̐ݒ�----------
		fout << "\nCELL_TYPES\t" << (this->inum - 1)*(this->jnum - 1) << "\n";
		for (int i = 0; i < (this->inum - 1)*(this->jnum - 1); i++) {
			fout << "9\n";
		}

		fout.close();
	}


	template<class T>
	void EllipticalMesher<T>::ExportForPANSFEM(std::string _fnameheader){
		//----------Node�̏o��----------
		std::ofstream fout_node(_fnameheader + "_Node.csv");

		//.....Header�̏o��.....
		fout_node << "ID,DOX,x0,x1";

		//.....���W���o��.....
		for (int i = 0; i < this->inum; i++) {
			for (int j = 0; j < this->jnum; j++) {
				fout_node << std::endl << i * (this->jnum) + j << "," << "2" << "," << this->pin[i][j].x[0] << "," << this->pin[i][j].x[1];
			}
		}

		fout_node.close();

		//----------Element�̏o��----------
		std::ofstream fout_element(_fnameheader + "_Element.csv");

		//.....Header�̏o��.....
		fout_element << "ID,NON,n0,n1,n2,n3";

		//.....�ߓ_�ԍ����o��.....
		for (int i = 0; i < this->inum - 1; i++) {
			for (int j = 0; j < this->jnum - 1; j++) {
				fout_element << std::endl << i * (this->jnum - 1) + j << "," << "4" << "," << i * (this->jnum) + j << "," << (i + 1)*(this->jnum) + j << "," << (i + 1)*(this->jnum) + j + 1 << "," << i * (this->jnum) + j + 1;
			}
		}

		fout_element.close();
	}
}