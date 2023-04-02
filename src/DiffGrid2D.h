#pragma once

#include "omp.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;

class DiffGrid2D
{
    double hx;
    double hy;
    MatrixXd U, Ux, Uy, Uxx, Uyy;
public:
    DiffGrid2D();
    DiffGrid2D(uint16_t Nx, uint16_t Ny, double hx, double hy);
    DiffGrid2D(const MatrixXd& U, double hx, double hy);
    DiffGrid2D(const DiffGrid2D& set_DiffGrid2D);
    ~DiffGrid2D();

    // value

    void SetValue(const MatrixXd& set_U);

    void Seth(const double set_hx, const double set_hy);

    const MatrixXd& GetValue();

    void AddOneEntry(uint16_t coor_x, uint16_t coor_y, double value) {
        U(coor_x, coor_y) += value;
    }

    // differential opeartor
    const MatrixXd& dx_f();

    const MatrixXd& dx_b();

    const MatrixXd& dy_f();

    const MatrixXd& dy_b();

    const MatrixXd& dxx();

    const MatrixXd& dyy();

    void PrintShape();
};

