#include "Source.h"
#include <iostream>
#include <vector>

using Eigen::VectorXd;
using Eigen::MatrixXd;

// source function
const double PI = atan(1.0)*4;

VectorXd Gaussian1D(int Nx, double dx, double sigma, double center) {
    VectorXd res(Nx);
    for (int i = 0; i < Nx; i++)
    {
        double xi = dx * i;
        res(i) = exp(-1 * (xi - center)*(xi - center) / 2 / (sigma*sigma));
    }
    return res;
}

VectorXd Sine1D(int Nx, double dx, double omega) {
    VectorXd res(Nx);
    for (int i = 0; i < Nx; i++)
    {
        double xi = dx * i;
        res(i) = sin(2* PI * omega * xi);
    }
    return res;
}

MatrixXd Gaussian2D(int Nx, int Ny, double dx, double dy, double sigma_x, double sigma_y, double center_x, double center_y) {
    MatrixXd G(Nx, Ny);
    double xi;
    double yi;
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            xi = dx * i;
            yi = dy * j;
            G(i,j) = exp( -1 * ((xi-center_x)*(xi-center_x)/(2*sigma_x*sigma_x) + (yi-center_y)*(yi-center_y)/(2*sigma_y*sigma_y) ) );
        }   
    }
    return G;
}

// Coor2D structure
std::ostream& operator<< (std::ostream& output, const Coor2D& c) {
    output << c.x << " " << c.y;
    return output;
}

// constructor
PointSource2D::PointSource2D() : N(0) {
    source_coor.reserve(N);
    source_fn.reserve(N);
}

PointSource2D::PointSource2D(uint16_t coor_x, uint16_t coor_y): N(1) {
    source_coor.push_back(Coor2D(coor_x, coor_y));
    source_fn.push_back(VectorXd::Zero(1));
}

PointSource2D::PointSource2D(uint16_t coor_x, uint16_t coor_y, const VectorXd& fn): N(1) {
    source_coor.push_back(Coor2D(coor_x, coor_y));
    source_fn.push_back(fn);
}

PointSource2D::PointSource2D(const std::vector<Coor2D>& coor) : N(coor.size()) {
    source_coor.reserve(N);
    source_fn.reserve(N);
    for (int i = 0; i < coor.size(); i++)
    {
        source_coor.emplace_back(coor[i]);
        source_fn.emplace_back(VectorXd::Zero(1));
    }
}

PointSource2D::PointSource2D(const std::vector<Coor2D>& coor, const VectorXd& fn) : N(coor.size()) {
    source_coor.reserve(N);
    source_fn.reserve(N);
    for (int i = 0; i < coor.size(); i++)
    {
        source_coor.emplace_back(coor[i]);
        source_fn.emplace_back(fn);
    }
}

PointSource2D::PointSource2D(const std::vector<Coor2D>& coor, const std::vector<VectorXd>& fn) : N(coor.size()) {
        source_coor.reserve(N);
        source_fn.reserve(N);
        for (int i = 0; i < coor.size(); i++)
        {
            source_coor.emplace_back(coor[i]);
            source_fn.emplace_back(fn[i]);
        }
    }

// member function

void PointSource2D::UpdatePML(uint16_t N_pml) {
    for (int i = 0; i < source_coor.size(); i++)
    {
        source_coor[i] += N_pml;
    }
}

uint16_t PointSource2D::GetSourceNum() {
    return N;
}

Coor2D PointSource2D::GetCoor(uint16_t source_idx) {
    return source_coor[source_idx];
}
VectorXd PointSource2D::GetFn(uint16_t source_idx) {
    return source_fn[source_idx];
}
double PointSource2D::GetValue(uint16_t source_idx, uint16_t time_idx) {
    return source_fn[source_idx][time_idx];
}

void PointSource2D::PrintInfo() {
    std::cout << "Source number: " << N << std::endl;
    for (int i = 0; i < N; i++)
    {
        std::cout << source_coor[i] << std::endl;
        std::cout << source_fn[i] << "\n" << std::endl;
    }
}