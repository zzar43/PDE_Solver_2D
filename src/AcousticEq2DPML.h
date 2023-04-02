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

class AcousticEq2DPML : public ModelAcousticEq2D
{
    uint16_t N_PML;
    double coef_PML;
    DiffGrid2D Vx0, Vx1, Vy0, Vy1, P0, P1, Q0, Q1, R0, R1;
    bool Record;
    std::vector<MatrixXd> SolData;

public:
    AcousticEq2DPML();
    ~AcousticEq2DPML() {}
    void SetRecord(bool r);
    void TimeUpdate(uint16_t idx_t);
    void Solve();
    void SaveSol(std::string savename);
    void SaveSolData(std::string savename, uint16_t step_t);
    MatrixXd GetSol() const;

};