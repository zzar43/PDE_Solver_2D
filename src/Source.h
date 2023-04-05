#pragma once

#include <iostream>
#include <cmath>
#include "Eigen/Dense"

using Eigen::VectorXf;
using Eigen::MatrixXf;

// source function
VectorXf Gaussian1D(int Nx, double dx, double sigma, double center);
VectorXf Sine1D(int Nx, double dx, double omega);
MatrixXf Gaussian2D(int Nx, int Ny, double dx, double dy, double sigma_x, double sigma_y, double center_x, double center_y);

struct Coor2D
{
    int x, y;
    Coor2D() : x(0), y(0) {}
    Coor2D(int x, int y) : x(x), y(y) {}
    int X() {
        return x;
    }
    int Y() {
        return y;
    }
    const Coor2D& operator+=(int N) {
        x += N;
        y += N;
        return *this;
    }
};

class PointSource2D
{
    // N: the number of point source
    int N = 0;
    std::vector<Coor2D> source_coor;
    std::vector<VectorXf> source_fn;
public:
    PointSource2D();
    PointSource2D(int coor_x, int coor_y);
    PointSource2D(int coor_x, int coor_y, const VectorXf& fn);
    PointSource2D(const std::vector<Coor2D>& coor);
    PointSource2D(const std::vector<Coor2D>& coor, const VectorXf& fn);
    PointSource2D(const std::vector<Coor2D>& coor, const std::vector<VectorXf>& fn);

    void UpdatePML(int N_pml);

    int GetSourceNum();
    Coor2D GetCoor(int source_idx);
    VectorXf GetFn(int source_idx);
    double GetValue(int source_idx, int time_idx);
    void PrintInfo();
};


