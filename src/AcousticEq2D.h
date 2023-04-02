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

class AcousticEq2D : public Model2D
{

    bool Record;
    MatrixXd A, B, C, Rho, Vx_init, Vy_init, P_init;
    PointSource2D source;
    DiffGrid2D Vx0, Vx1, Vy0, Vy1, P0, P1;

public:
    AcousticEq2D();
    ~AcousticEq2D() {}

    void SetParameters(const MatrixXd& set_C, const MatrixXd& set_Rho);
    void SetInit(const MatrixXd& set_Vx_init, const MatrixXd& set_Vy_init, const MatrixXd& set_P_init);
    void SetSource(const PointSource2D& set_source);
    void SetRecord(bool r);

    void TimeUpdate(uint16_t idx_t);
    void Solve();

    void SaveSol(std::string savename);
    void SaveSolData(std::string savename, uint16_t step_t);

    MatrixXd GetSol() const;
};