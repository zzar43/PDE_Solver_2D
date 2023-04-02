#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "omp.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;

class DiffGrid2D
{
    double hx;
    double hy;
    MatrixXd U, Ux, Uy, Uxx, Uyy;
public:
    DiffGrid2D();
    DiffGrid2D(uint16_t Nx, uint16_t Ny, double hx, double hy);
    DiffGrid2D(const MatrixXd& U, double hx, double hy);
    DiffGrid2D(const DiffGrid2D& set_DiffGrid2D);
    ~DiffGrid2D();

    // value

    void SetValue(const MatrixXd& set_U);

    void Seth(const double set_hx, const double set_hy);

    const MatrixXd& GetValue();

    void AddOneEntry(uint16_t coor_x, uint16_t coor_y, double value) {
        U(coor_x, coor_y) += value;
    }

    // differential opeartor
    const MatrixXd& dx_f();
    const MatrixXd& dx_b();
    const MatrixXd& dy_f();
    const MatrixXd& dy_b();
    const MatrixXd& dxx();
    const MatrixXd& dyy();

    void PrintShape();
};

// a general 2D model class for defining basic parameters
struct Model2D
{
    uint16_t Nx, Ny, N_pml, Nt;
    double hx, hy, tau, coef_pml;
    MatrixXd sigmaX, sigmaY;
    // save result
    MatrixXd Res;
    std::vector<MatrixXd> AllRes;

    // constructor
    Model2D();
    
    // member functions
    // set values
    void SetSpatial(uint16_t set_Nx, uint16_t set_Ny, const double& set_hx, const double& set_hy);
    void SetTime(uint16_t set_Nt, const double& set_tau);
    void SetPML(uint16_t set_N_pml, const double& set_coef_pml);

    // PML functions
    void BuildSigmaPML();
    MatrixXd ExpandToPML(const MatrixXd& M);
    MatrixXd ContractFromPML(const MatrixXd& M_pml);

    // Print
    void PrintInfo();

    // save
    void SaveRes(std::string savename);
    void SaveAllRes(std::string savename, uint16_t step_t);

};