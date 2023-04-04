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

using Eigen::MatrixXd;

class MaxwellsEqTE2DPML : public Model2D
{
    bool Record;
    MatrixXd A, B, Eps, Mu, Ex_init, Ey_init, Hz_init;
    PointSource2D source;
    DiffGrid2D Ex0, Ex1, Ey0, Ey1, Hz0, Hz1, P0, P1, Q0, Q1, R0, R1;

public:
    MaxwellsEqTE2DPML();
    ~MaxwellsEqTE2DPML() {}

    void SetParameters(const MatrixXd& set_Eps, const MatrixXd& set_Mu);
    void SetInit(const MatrixXd& set_Ex_init, const MatrixXd& set_Ey_init, const MatrixXd& set_Hz_init);
    void SetSource(const PointSource2D& set_source);
    void SetRecord(bool r);

    void TimeUpdate(uint16_t idx_t);
    void Solve();

    MatrixXd GetSol() const;

};