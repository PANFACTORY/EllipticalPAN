//*****************************************************************************
// Title		:src/cpp/EllipticalMesher.h
// Author	:Tanabe Yuta
// Date		:2019/10/25
// Copyright	:(C)2019 TanabeYuta
//*****************************************************************************

#pragma once
#include <iostream>

namespace EllipticalPAN {
template <template <class> class P, class T>
class Mesher {
   public:
    Mesher() = delete;
    Mesher(int _ni, int _nj, T _DEPS = T(1e-8), int _ITRMAX = 100000)
        : ni(_ni), nj(_nj), DEPS(_DEPS), ITRMAX(_ITRMAX) {
        this->p = new P<T>[this->ni * this->nj];
    }
    Mesher(const Mesher<P, T>&) = delete;
    ~Mesher() { delete[] this->p; }

    void SetPoint(int _i, int _j, const P<T>& point) {
        int i = _i == -1 ? this->ni - 1 : _i, j = _j == -1 ? this->nj - 1 : _j;
        this->p[this->ID(i, j)] = point;
    }
    int Generate() {
        //----------Solve Laplace equation----------
        for (int itr = 0; itr < this->ITRMAX; itr++) {
            T dxmax = T();
            for (int i = 1; i < this->ni - 1; i++) {
                for (int j = 1; j < this->nj - 1; j++) {
                    //.....Make Laplace equation.....
                    P<T> xxi = 0.5 * (this->p[this->ID(i + 1, j)] -
                                      this->p[this->ID(i - 1, j)]);
                    P<T> xiota = 0.5 * (this->p[this->ID(i, j + 1)] -
                                        this->p[this->ID(i, j - 1)]);
                    T alpha = xiota.Dot(xiota) + 1.0e-8;
                    T beta = xxi.Dot(xiota);
                    T ganma = xxi.Dot(xxi) + 1.0e-8;

                    //.....Update values.....
                    P<T> tmp = this->p[this->ID(i, j)];
                    this->p[this->ID(i, j)] =
                        (0.5 / (alpha + ganma)) *
                        (alpha * (this->p[this->ID(i + 1, j)] +
                                  this->p[this->ID(i - 1, j)]) -
                         0.5 * beta *
                             (this->p[this->ID(i + 1, j + 1)] -
                              this->p[this->ID(i - 1, j + 1)] -
                              this->p[this->ID(i + 1, j - 1)] +
                              this->p[this->ID(i - 1, j - 1)]) +
                         ganma * (this->p[this->ID(i, j + 1)] +
                                  this->p[this->ID(i, j - 1)]));

                    //.....Get norm maximam.....
                    T dxtmp = (tmp - this->p[this->ID(i, j)]).Norm();
                    if (dxmax < dxtmp) {
                        dxmax = dxtmp;
                    }
                }
            }
            //.....Check convergence.....
            if (dxmax < DEPS) {
                return itr + 1;
            }
        }
        return this->ITRMAX;
    }
    void ExportAsVTK(std::ostream& fout) {
        //----------Headers----------
        fout << "# vtk DataFile Version 4.1" << std::endl;
        fout << "vtk output" << std::endl;
        fout << "ASCII" << std::endl;
        fout << "DATASET UNSTRUCTURED_GRID" << std::endl;

        //----------Points----------
        fout << std::endl
             << "POINTS\t" << this->ni * this->nj << "\tfloat" << std::endl;
        for (int i = 0; i < this->ni; i++) {
            for (int j = 0; j < this->nj; j++) {
                fout << this->p[this->ID(i, j)][0] << "\t"
                     << this->p[this->ID(i, j)][1] << "\t" << 0.0 << std::endl;
            }
        }

        //----------Cells----------
        fout << std::endl
             << "CELLS " << (this->ni - 1) * (this->nj - 1) << "\t"
             << (this->ni - 1) * (this->nj - 1) * 5 << std::endl;
        for (int i = 0; i < this->ni - 1; i++) {
            for (int j = 0; j < this->nj - 1; j++) {
                fout << "4 " << i * (this->nj) + j << "\t"
                     << i * (this->nj) + j + 1 << "\t"
                     << (i + 1) * (this->nj) + j + 1 << "\t"
                     << (i + 1) * (this->nj) + j << std::endl;
            }
        }

        //----------Cell types----------
        fout << std::endl
             << "CELL_TYPES\t" << (this->ni - 1) * (this->nj - 1) << std::endl;
        for (int i = 0; i < (this->ni - 1) * (this->nj - 1); i++) {
            fout << "9" << std::endl;
        }
    }

   private:
    const int ni, nj, ITRMAX;
    const T DEPS;

    P<T>* p;

    int ID(int i, int j) const { return i + this->ni * j; }
};
}  // namespace EllipticalPAN
