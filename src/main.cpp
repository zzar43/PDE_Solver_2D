#include <iostream>
#include <chrono>

#include "Eigen/Dense"
#include "DiffGrid2D.h"
// #include "WaveEq2D.h"
#include "AcousticEq2DPML.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using namespace chrono; 

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

int main() {

    AcousticEq2DPML model;

    model.SetSpatial(201,201,0.005,0.005);
    model.SetTime(1001,0.0005);
    model.SetPML(30,0.1);
    
    MatrixXd cc = MatrixXd::Ones(model.Nx, model.Ny);
    MatrixXd rho = MatrixXd::Ones(model.Nx, model.Ny);
    MatrixXd vx = MatrixXd::Zero(model.Nx, model.Ny);
    MatrixXd vy = MatrixXd::Zero(model.Nx, model.Ny);
    MatrixXd pp = MatrixXd::Zero(model.Nx, model.Ny);
    // MatrixXd pp = Gaussian2D(model.Nx, model.Ny, model.hx, model.hy, 0.1, 0.1, 0.5, 0.5);
    model.SetParameters(cc, rho);
    model.SetInit(vx, vy, pp);

    // source
    Coor2D c1(51,51);
    Coor2D c2(51,151);
    Coor2D c3(151,51);
    Coor2D c4(151,151);
    vector<Coor2D> c{c1,c2,c3,c4};
    VectorXd f1 = Sine1D(model.Nt, model.tau, 5);
    VectorXd f2 = Sine1D(model.Nt, model.tau, 10);
    VectorXd f3 = Sine1D(model.Nt, model.tau, 15);
    VectorXd f4 = Sine1D(model.Nt, model.tau, 20);
    vector<VectorXd> f{f1,f2,f3,f4};
    PointSource2D ptsource(c, f);
    model.SetSource(ptsource);

    // solve
    model.SetRecord(true);

    auto time1 = high_resolution_clock::now();
    model.Solve();
    auto time2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(time2 - time1);
    cout << "Computation time: " << duration.count() << " milliseconds"<< endl;

    // record
    model.SaveRes(".\\data\\sol.txt");
    model.SaveAllRes(".\\data\\solData.txt", 100);

}