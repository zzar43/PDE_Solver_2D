#include "Source.h"
#include <iostream>
#include <vector>

using Eigen::VectorXf;
using Eigen::MatrixXf;

// source function
const double PI = atan(1.0)*4;

VectorXf Gaussian1D(int Nx, double dx, double sigma, double center) {
    VectorXf res(Nx);
    for (int i = 0; i < Nx; i++)
    {
        double xi = dx * i;
        res(i) = exp(-1 * (xi - center)*(xi - center) / 2 / (sigma*sigma));
    }
    return res;
}

VectorXf Sine1D(int Nx, double dx, double omega) {
    VectorXf res(Nx);
    for (int i = 0; i < Nx; i++)
    {
        double xi = dx * i;
        res(i) = sin(2* PI * omega * xi);
    }
    return res;
}

MatrixXf Gaussian2D(int Nx, int Ny, double dx, double dy, double sigma_x, double sigma_y, double center_x, double center_y) {
    MatrixXf G(Nx, Ny);
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

PointSource2D::PointSource2D(int coor_x, int coor_y): N(1) {
    source_coor.push_back(Coor2D(coor_x, coor_y));
    source_fn.push_back(VectorXf::Zero(1));
}

PointSource2D::PointSource2D(int coor_x, int coor_y, const VectorXf& fn): N(1) {
    source_coor.push_back(Coor2D(coor_x, coor_y));
    source_fn.push_back(fn);
}

PointSource2D::PointSource2D(const std::vector<Coor2D>& coor) : N(coor.size()) {
    source_coor.reserve(N);
    source_fn.reserve(N);
    for (int i = 0; i < coor.size(); i++)
    {
        source_coor.emplace_back(coor[i]);
        source_fn.emplace_back(VectorXf::Zero(1));
    }
}

PointSource2D::PointSource2D(const std::vector<Coor2D>& coor, const VectorXf& fn) : N(coor.size()) {
    source_coor.reserve(N);
    source_fn.reserve(N);
    for (int i = 0; i < coor.size(); i++)
    {
        source_coor.emplace_back(coor[i]);
        source_fn.emplace_back(fn);
    }
}

PointSource2D::PointSource2D(const std::vector<Coor2D>& coor, const std::vector<VectorXf>& fn) : N(coor.size()) {
        source_coor.reserve(N);
        source_fn.reserve(N);
        for (int i = 0; i < coor.size(); i++)
        {
            source_coor.emplace_back(coor[i]);
            source_fn.emplace_back(fn[i]);
        }
    }

// member function

void PointSource2D::UpdatePML(int N_pml) {
    for (int i = 0; i < source_coor.size(); i++)
    {
        source_coor[i] += N_pml;
    }
}

int PointSource2D::GetSourceNum() {
    return N;
}

Coor2D PointSource2D::GetCoor(int source_idx) {
    return source_coor[source_idx];
}
VectorXf PointSource2D::GetFn(int source_idx) {
    return source_fn[source_idx];
}
double PointSource2D::GetValue(int source_idx, int time_idx) {
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