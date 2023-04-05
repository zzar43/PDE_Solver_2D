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
    uint16_t x, y;
    Coor2D() : x(0), y(0) {}
    Coor2D(uint16_t x, uint16_t y) : x(x), y(y) {}
    uint16_t X() {
        return x;
    }
    uint16_t Y() {
        return y;
    }
    const Coor2D& operator+=(uint16_t N) {
        x += N;
        y += N;
        return *this;
    }
};

class PointSource2D
{
    // N: the number of point source
    uint16_t N = 0;
    std::vector<Coor2D> source_coor;
    std::vector<VectorXf> source_fn;
public:
    PointSource2D();
    PointSource2D(uint16_t coor_x, uint16_t coor_y);
    PointSource2D(uint16_t coor_x, uint16_t coor_y, const VectorXf& fn);
    PointSource2D(const std::vector<Coor2D>& coor);
    PointSource2D(const std::vector<Coor2D>& coor, const VectorXf& fn);
    PointSource2D(const std::vector<Coor2D>& coor, const std::vector<VectorXf>& fn);

    void UpdatePML(uint16_t N_pml);

    uint16_t GetSourceNum();
    Coor2D GetCoor(uint16_t source_idx);
    VectorXf GetFn(uint16_t source_idx);
    double GetValue(uint16_t source_idx, uint16_t time_idx);
    void PrintInfo();
};


