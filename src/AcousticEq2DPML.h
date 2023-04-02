#pragma once

// Consider the acoustic wave equation in the first order form in a 2D plane with PML.

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Eigen/Dense"
#include "DiffGrid2D.h"
#include "Source.h"
#include "AcousticEq2D.h"

using Eigen::MatrixXd;

class AcousticEq2DPML : public Model2D
{
    bool Record;
    MatrixXd A, B, C, Rho, Vx_init, Vy_init, P_init;
    PointSource2D source;
    DiffGrid2D Vx0, Vx1, Vy0, Vy1, P0, P1, Q0, Q1, R0, R1;

public:
    AcousticEq2DPML();
    ~AcousticEq2DPML() {}

    void SetParameters(const MatrixXd& set_C, const MatrixXd& set_Rho);
    void SetInit(const MatrixXd& set_Vx_init, const MatrixXd& set_Vy_init, const MatrixXd& set_P_init);
    void SetSource(const PointSource2D& set_source);
    void SetRecord(bool r);

    void TimeUpdate(uint16_t idx_t);
    void Solve();

    MatrixXd GetSol() const;

};