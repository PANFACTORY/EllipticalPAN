//*****************************************************************************
// Title		:src/cpp/EllipticalMesher.h
// Author	:Tanabe Yuta
// Date		:2019/10/25
// Copyright	:(C)2019 TanabeYuta
//*****************************************************************************

#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "point.h"

namespace EllipticalPAN {
template <class T>
class Mesher {
   public:
    Mesher() = delete;
    Mesher(int _inum, int _jnum, T _DEPS = T(1e-8), int _ITRMAX = 100000)
        : inum(_inum), jnum(_jnum), DEPS(_DEPS), ITRMAX(_ITRMAX) {
        this->pin = std::vector<std::vector<Point<T> > >(
            this->inum, std::vector<Point<T> >(this->jnum, Point<T>(T(), T())));
    }
    Mesher(const Mesher<T>&) = delete;
    ~Mesher() {}

    void SetPoint(int _i, int _j, const Point<T>& point) {
        int i = _i == -1 ? this->inum - 1 : _i,
            j = _j == -1 ? this->jnum - 1 : _j;
        this->pin[i][j] = point;
    }
    void Generate() {
        //----------Solve Laplace equation----------
        for (int k = 0; k < this->ITRMAX; k++) {
            T dxmax = T();
            for (int i = 1; i < this->inum - 1; i++) {
                for (int j = 1; j < this->jnum - 1; j++) {
                    //.....Make Laplace equation.....
                    Point<T> xxi =
                        0.5 * (this->pin[i + 1][j] - this->pin[i - 1][j]);
                    Point<T> xiota =
                        0.5 * (this->pin[i][j + 1] - this->pin[i][j - 1]);
                    T alpha = xiota.Dot(xiota) + 1.0e-8;
                    T beta = xxi.Dot(xiota);
                    T ganma = xxi.Dot(xxi) + 1.0e-8;

                    //.....Update values.....
                    Point<T> tmp = this->pin[i][j];
                    this->pin[i][j] =
                        (0.5 / (alpha + ganma)) *
                        (alpha * (this->pin[i + 1][j] + this->pin[i - 1][j]) -
                         0.5 * beta *
                             (this->pin[i + 1][j + 1] -
                              this->pin[i - 1][j + 1] -
                              this->pin[i + 1][j - 1] +
                              this->pin[i - 1][j - 1]) +
                         ganma * (this->pin[i][j + 1] + this->pin[i][j - 1]));

                    //.....Get norm maximam.....
                    T dxtmp = (tmp - this->pin[i][j]).Norm();
                    if (dxmax < dxtmp) {
                        dxmax = dxtmp;
                    }
                }
            }
            //.....Check convergence.....
            if (dxmax < DEPS) {
                std::cout << "Converged:k = " << k + 1;
                break;
            }
        }
    }
    void ExportToVTK(std::string _fname) {
        std::ofstream fout(_fname + ".vtk");

        //----------Headers----------
        fout << "# vtk DataFile Version 4.1\n";
        fout << "vtk output\n";
        fout << "ASCII\n";
        fout << "DATASET UNSTRUCTURED_GRID\n";

        //----------Points----------
        fout << "\nPOINTS\t" << this->inum * this->jnum << "\tfloat\n";
        for (auto pinrow : this->pin) {
            for (auto pinone : pinrow) {
                fout << pinone[0] << "\t" << pinone[1] << "\t" << 0.0 << "\n";
            }
        }

        //----------Cells----------
        fout << "\nCELLS " << (this->inum - 1) * (this->jnum - 1) << "\t"
             << (this->inum - 1) * (this->jnum - 1) * 5 << "\n";
        for (int i = 0; i < this->inum - 1; i++) {
            for (int j = 0; j < this->jnum - 1; j++) {
                fout << "4 " << i * (this->jnum) + j << "\t"
                     << i * (this->jnum) + j + 1 << "\t"
                     << (i + 1) * (this->jnum) + j + 1 << "\t"
                     << (i + 1) * (this->jnum) + j << "\n";
            }
        }

        //----------Cell types----------
        fout << "\nCELL_TYPES\t" << (this->inum - 1) * (this->jnum - 1) << "\n";
        for (int i = 0; i < (this->inum - 1) * (this->jnum - 1); i++) {
            fout << "9\n";
        }

        fout.close();
    }

   private:
    const int inum, jnum;

    std::vector<std::vector<Point<T> > > pin;

    const int ITRMAX;
    const T DEPS;
};
}  // namespace EllipticalPAN
