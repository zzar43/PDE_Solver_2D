#pragma once

// Classic wave equation in 2D spatial space: \frac{1}{c^2} u_{tt} - \Delta u = f
// Finite difference scheme, central difference in space

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Eigen/Dense"
#include "DiffGrid2D.h"
#include "Source.h"

using Eigen::MatrixXd;

struct ModelWaveEq2D
{
    uint16_t Nx, Ny, Nt;
    double hx, hy, tau;
    MatrixXd C, C2, U_init;
    PointSource2D source;

    // constructor
    ModelWaveEq2D();
    ModelWaveEq2D(uint16_t Nx, uint16_t Ny, uint16_t Nt, double hx, double hy, double tau, MatrixXd C, MatrixXd U_init)
        : Nx(Nx), Ny(Ny), Nt(Nt), hx(hx), hy(hy), tau(tau), C(C), U_init(U_init) {}

    // set model value
    void SetSpatial(uint16_t set_Nx, uint16_t set_Ny, const double& set_hx, const double& set_hy);
    void SetTime(uint16_t set_Nt, const double& set_tau);
    void SetC(const MatrixXd& set_C);
    void SetInit(const MatrixXd& set_U_init);
    void SetSource(const PointSource2D& set_source);
    // no set boundary now, just use Dirichlet boundary condition

    void PrintInfo() const;
};


class WaveEq2D : public ModelWaveEq2D
{
    MatrixXd LU;
    DiffGrid2D U0, U1, U2;
    bool Record;
    std::vector<MatrixXd> SolData;

public:
    WaveEq2D();
    ~WaveEq2D() {}
    void SetRecord(bool r);
    void EvalLU(DiffGrid2D& U);
    void TimeUpdate(uint16_t idx_t);
    void Solve();
    void SaveSol(std::string savename);
    void SaveSolData(std::string savename, uint16_t step_t);
    MatrixXd GetSol();
    MatrixXd GetLU();
};