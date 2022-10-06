/**
 * @file mesher.h
 * @author PANFACTORY (github/PANFACTORY)
 * @brief Define and implement the Mesher class
 * @version 0.1
 * @date 2022-10-05
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#include <iostream>

namespace EllipticalPAN {
/**
 * @brief Mesher class
 *
 * @tparam P    Template class of Point
 *              This must contain thise methods:
 *                  Default constructor setting x and y to 0
 *                  Copy constructor
 *                  operator+ : vector vector addition
 *                  operator- : vector vector subtruction
 *                  operator* : vector scalar multiplication
 *                  operator/ : vector scalar division
 *                  operator(): access each element for example (0) is x
 *                  dot(P)    : scalar product with P
 *                  norm()    : norm of the vector
 * @tparam T    Type of variable for example double
 */
template <class P, class T>
class Mesher {
   public:
    Mesher() = delete;
    Mesher(const Mesher<P, T>&) = delete;
    ~Mesher() { delete[] this->p; }

    /**
     * @brief Construct a new Mesher object
     *
     * @param _ni       Number of horizontal grid points
     * @param _nj       Number of vertical grid points
     */
    Mesher(int _ni, int _nj) : ni(_ni), nj(_nj) {
        this->p = new P[this->ni * this->nj];
    }

    /**
     * @brief Set the Point object
     *
     * @param _i    Horizontal position to set
     * @param _j    Vertical position to set
     * @param point Point object set
     */
    void SetPoint(int _i, int _j, const P& point) {
        int i = _i == -1 ? this->ni - 1 : _i, j = _j == -1 ? this->nj - 1 : _j;
        this->p[this->ID(i, j)] = point;
    }

    /**
     * @brief Get the Point object
     *
     * @param _i    Horizontal position to get
     * @param _j    Vertical position to get
     * @return P    Point object to get
     */
    P GetPoint(int _i, int _j) const {
        int i = _i == -1 ? this->ni - 1 : _i, j = _j == -1 ? this->nj - 1 : _j;
        return this->p[this->ID(i, j)];
    }

    /**
     * @brief Generate mesh
     *
     * @param DEPS      Convergence criterion
     * @param ITRMAX    Maximum number of iteration
     * @return int      Number of steps at convergence
     */
    int Generate(T DEPS = T(1e-8), int ITRMAX = 100000) {
        //----------Solve Laplace equation----------
        for (int itr = 0; itr < ITRMAX; itr++) {
            T dxmax = T();
            for (int i = 1; i < this->ni - 1; i++) {
                for (int j = 1; j < this->nj - 1; j++) {
                    //.....Make Laplace equation.....
                    P xxi = 0.5 * (this->p[this->ID(i + 1, j)] -
                                   this->p[this->ID(i - 1, j)]);
                    P xiota = 0.5 * (this->p[this->ID(i, j + 1)] -
                                     this->p[this->ID(i, j - 1)]);
                    T alpha = xiota.dot(xiota) + T(1.0e-8);
                    T beta = xxi.dot(xiota) + T(1.0e-8);
                    T ganma = xxi.dot(xxi) + T(1.0e-8);

                    //.....Update values.....
                    P tmp = this->p[this->ID(i, j)];
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
                    T dxtmp = (tmp - this->p[this->ID(i, j)]).norm();
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
        return ITRMAX;
    }

    /**
     * @brief Export mesh as vtk format to stream
     *
     * @param fout  Stream updated
     */
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
                fout << this->p[this->ID(i, j)](0) << "\t"
                     << this->p[this->ID(i, j)](1) << "\t" << T() << std::endl;
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
    const int ni, nj;

    P* p;

    int ID(int i, int j) const { return i + this->ni * j; }
};
}  // namespace EllipticalPAN
