#pragma once

// Consider the convection-diffusion equation in 2D
// u_t = \nabla \cdot (D\nabla u) - \nabla \cdot (vu) + R

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Eigen/Dense"
#include "DiffGrid2D.h"
#include "Source.h"

using Eigen::MatrixXf;

class ConvDiffEq2D : public Model2D
{
    bool Record;
    MatrixXf D, Vx, Vy, U_init;
    PointSource2D source;
    DiffGrid2D U0, U1, DdxU, DdyU, VxU, VyU;

public:
    ConvDiffEq2D();
    ~ConvDiffEq2D() {}

    void SetParameters(const MatrixXf& set_D, const MatrixXf set_Vx, const MatrixXf set_Vy);
    void SetInit(const MatrixXf& set_U_init);
    void SetSource(const PointSource2D& set_source);
    void SetRecord(bool r);

    void TimeUpdate(int idx_t);
    void Solve();

    MatrixXf GetSol() const;
};