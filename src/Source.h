#pragma once

#include <iostream>
#include <cmath>
#include "Eigen/Dense"

using Eigen::VectorXd;

// source function
VectorXd Gaussian1D(int Nx, double dx, double sigma, double center);
VectorXd Sine1D(int Nx, double dx, double omega);

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
};

class PointSource2D
{
    // N: the number of point source
    uint16_t N = 0;
    std::vector<Coor2D> source_coor;
    std::vector<VectorXd> source_fn;
public:
    PointSource2D();
    PointSource2D(uint16_t coor_x, uint16_t coor_y);
    PointSource2D(uint16_t coor_x, uint16_t coor_y, const VectorXd& fn);
    PointSource2D(const std::vector<Coor2D>& coor);
    PointSource2D(const std::vector<Coor2D>& coor, const VectorXd& fn);
    PointSource2D(const std::vector<Coor2D>& coor, const std::vector<VectorXd>& fn);

    uint16_t GetSourceNum();
    Coor2D GetCoor(uint16_t source_idx);
    VectorXd GetFn(uint16_t source_idx);
    double GetValue(uint16_t source_idx, uint16_t time_idx);
    void PrintInfo();
};


