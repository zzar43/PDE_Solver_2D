#pragma once

// Consider the Maxwell's equations in a lossless and homogeneous medium.
// Consider the transverse electric (TE) mode.
// \eps \partial_t E_x = \partial_y H_z
// \eps \partial_t E_y = - \partial_x H_z
// \mu \partial_t H_z = \partial_y E_x - \partial_x E_y

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Eigen/Dense"
#include "DiffGrid2D.h"
#include "Source.h"

using Eigen::MatrixXf;

class MaxwellsEqTE2DPML : public Model2D
{
    bool Record;
    MatrixXf A, B, Eps, Mu, Ex_init, Ey_init, Hz_init;
    PointSource2D source;
    DiffGrid2D Ex0, Ex1, Ey0, Ey1, Hz0, Hz1, P0, P1, Q0, Q1, R0, R1;

public:
    MaxwellsEqTE2DPML();
    ~MaxwellsEqTE2DPML() {}

    void SetParameters(const MatrixXf& set_Eps, const MatrixXf& set_Mu);
    void SetInit(const MatrixXf& set_Ex_init, const MatrixXf& set_Ey_init, const MatrixXf& set_Hz_init);
    void SetSource(const PointSource2D& set_source);
    void SetRecord(bool r);

    void TimeUpdate(int idx_t);
    void Solve();

    MatrixXf GetSol() const;

};