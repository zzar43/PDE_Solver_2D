#include <iostream>
#include "DiffGrid2D.h"

using Eigen::MatrixXd;

DiffGrid2D::DiffGrid2D() {
    hx = hy = 0.;
    U = MatrixXd();
    Ux = MatrixXd();
    Uy = MatrixXd();
    Uxx = MatrixXd();
    Uyy = MatrixXd();
}

DiffGrid2D::DiffGrid2D(uint16_t Nx, uint16_t Ny, double hx, double hy) : hx(hx), hy(hy)
{
    U = MatrixXd::Zero(Nx, Ny);
    Ux = MatrixXd::Zero(Nx, Ny);
    Uy = MatrixXd::Zero(Nx, Ny);
    Uxx = MatrixXd::Zero(Nx, Ny);
    Uyy = MatrixXd::Zero(Nx, Ny);
}

DiffGrid2D::DiffGrid2D(const MatrixXd& set_U, double hx, double hy) : hx(hx), hy(hy)
{
    U = set_U;
    Ux = MatrixXd::Zero(U.rows(), U.cols());
    Uy = MatrixXd::Zero(U.rows(), U.cols());
    Uxx = MatrixXd::Zero(U.rows(), U.cols());
    Uyy = MatrixXd::Zero(U.rows(), U.cols());
}

DiffGrid2D::DiffGrid2D(const DiffGrid2D& set_DiffGrid2D) {
    hx = set_DiffGrid2D.hx;
    hy = set_DiffGrid2D.hy;
    U = set_DiffGrid2D.U;
    Ux = set_DiffGrid2D.Ux;
    Uy = set_DiffGrid2D.Uy;
    Uxx = set_DiffGrid2D.Uxx;
    Uyy = set_DiffGrid2D.Uyy;
}

DiffGrid2D::~DiffGrid2D() {}

void DiffGrid2D::SetValue(const MatrixXd& set_U) {
    U = set_U;
}

void DiffGrid2D::Seth(const double set_hx, const double set_hy) {
    hx = set_hx;
    hy = set_hy;
}

const MatrixXd& DiffGrid2D::GetValue() {
    return U;
}

// differential opeartor, forward and backward
const MatrixXd& DiffGrid2D::dx_f() {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < U.cols()-1; j++)
    {
        for (int i = 1; i < U.rows()-1; i++)
        {
            Ux(i,j) = (U(i+1,j) - U(i,j)) / hx;
        }
    }
    return Ux;
}

const MatrixXd& DiffGrid2D::dx_b() {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < U.cols()-1; j++)
    {
        for (int i = 1; i < U.rows()-1; i++)
        {
            Ux(i,j) = (U(i,j) - U(i-1,j)) / hx;
        }
    }
    return Ux;
}

const MatrixXd& DiffGrid2D::dy_f() {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < U.cols()-1; j++)
    {
        for (int i = 1; i < U.rows()-1; i++)
        {
            Uy(i,j) = (U(i,j+1) - U(i,j)) / hy;
        }
    }
    return Uy;
}

const MatrixXd& DiffGrid2D::dy_b() {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < U.cols()-1; j++)
    {
        for (int i = 1; i < U.rows()-1; i++)
        {
            Uy(i,j) = (U(i,j) - U(i,j-1)) / hy;
        }
    }
    return Uy;
}

const MatrixXd& DiffGrid2D::dxx() {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < U.cols()-1; j++)
    {
        for (int i = 1; i < U.rows()-1; i++)
        {
            Uxx(i,j) = (U(i+1,j) - 2*U(i,j) + U(i-1,j)) / hx / hx;
        }
    }
    return Uxx;
}

const MatrixXd& DiffGrid2D::dyy() {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < U.cols()-1; j++)
    {
        for (int i = 1; i < U.rows()-1; i++)
        {
            Uyy(i,j) = (U(i,j+1) - 2*U(i,j) + U(i,j-1)) / hy / hy;
        }
    }
    return Uyy;
}

void DiffGrid2D::PrintShape() {
    std::cout << " Nx = " << U.rows() << " Ny = " << U.cols();
    std::cout << " hx = " << hx << " hy = " << hy << std::endl;
}
