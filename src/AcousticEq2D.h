#pragma once

// Consider the acoustic wave equation in the first order form:
// \frac{\partial v}{\partial t} = - 1 / rho_0 \nabla p
// \frac{\partial p}{\partial t} = - \kappa_0 \nabla\cdot v + f
// where \rho_0 is the mass density, \kappa_0 = \rho_0 c^2 is the bulk modulus, c is the wave speed.
// The source is add onto p
// Let A = - 1 / rho_0, B = - \kappa_0
// Finite difference scheme is used, simple upwind scheme is applied

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Eigen/Dense"
#include "DiffGrid2D.h"
#include "Source.h"

using Eigen::MatrixXd;

struct ModelAcousticEq2D
{
    uint16_t Nx, Ny, Nt;
    double hx, hy, tau;
    MatrixXd A, B, C, Rho, Vx_init, Vy_init, P_init;
    PointSource2D source;

    // constructor
    ModelAcousticEq2D();
    ModelAcousticEq2D(uint16_t Nx, uint16_t Ny, uint16_t Nt, double hx, double hy, double tau, MatrixXd C, MatrixXd Rho, MatrixXd Vx_init, MatrixXd Vy_init, MatrixXd P_init)
        : Nx(Nx), Ny(Ny), Nt(Nt), hx(hx), hy(hy), tau(tau), C(C), Rho(Rho), Vx_init(Vx_init), Vy_init(Vy_init), P_init(P_init) {}

    // set model value
    void SetSpatial(uint16_t set_Nx, uint16_t set_Ny, const double& set_hx, const double& set_hy);
    void SetTime(uint16_t set_Nt, const double& set_tau);
    void SetModel(const MatrixXd& set_C, const MatrixXd& set_Rho);
    void SetInit(const MatrixXd& set_Vx_init, const MatrixXd& set_Vy_init, const MatrixXd& set_P_init);
    void SetSource(const PointSource2D& set_source);
    // no set boundary now, just use Dirichlet boundary condition

    void PrintInfo() const;
};

class AcousticEq2D : public ModelAcousticEq2D
{
    DiffGrid2D Vx0, Vx1, Vy0, Vy1, P0, P1;
    bool Record;
    std::vector<MatrixXd> SolData;

public:
    AcousticEq2D();
    ~AcousticEq2D() {}
    void SetRecord(bool r);
    void TimeUpdate(uint16_t idx_t);
    void Solve();
    void SaveSol(std::string savename);
    void SaveSolData(std::string savename, uint16_t step_t);
    MatrixXd GetSol() const;
};