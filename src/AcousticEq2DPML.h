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

using Eigen::MatrixXf;

class AcousticEq2DPML : public Model2D
{
    bool Record;
    MatrixXf A, B, C, Rho, Vx_init, Vy_init, P_init;
    PointSource2D source;
    DiffGrid2D Vx0, Vx1, Vy0, Vy1, P0, P1, Q0, Q1, R0, R1;

public:
    AcousticEq2DPML();
    ~AcousticEq2DPML() {}

    void SetParameters(const MatrixXf& set_C, const MatrixXf& set_Rho);
    void SetInit(const MatrixXf& set_Vx_init, const MatrixXf& set_Vy_init, const MatrixXf& set_P_init);
    void SetSource(const PointSource2D& set_source);
    void SetRecord(bool r);

    void TimeUpdate(int idx_t);
    void Solve();

    MatrixXf GetSol() const;

};