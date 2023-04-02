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


class WaveEq2D : public Model2D
{
    bool Record;
    MatrixXd LU, C, C2, U_init;
    DiffGrid2D U0, U1, U2;
    PointSource2D source;
    // std::vector<MatrixXd> SolData;

public:
    WaveEq2D();
    ~WaveEq2D() {}
    void SetParameters(const MatrixXd& set_C);
    void SetInit(const MatrixXd& set_U_init);
    void SetSource(const PointSource2D& set_source);
    void SetRecord(bool r);
    void EvalLU(DiffGrid2D& U);
    void TimeUpdate(uint16_t idx_t);
    void Solve();
};