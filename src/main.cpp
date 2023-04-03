#include <iostream>
#include <chrono>

#include "Eigen/Dense"
#include "DiffGrid2D.h"
// #include "WaveEq2D.h"
// #include "AcousticEq2DPML.h"
#include "MaxwellsEqTE2DPML.h"

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

    MaxwellsEqTE2DPML model;

    model.SetSpatial(201,201,0.005,0.005);
    model.SetTime(301,0.002);
    model.SetPML(30,0.1);
    
    MatrixXd mu = MatrixXd::Ones(model.Nx, model.Ny);
    // MatrixXd eps = MatrixXd::Ones(model.Nx, model.Ny);
    MatrixXd eps = 1*Gaussian2D(model.Nx, model.Ny, model.hx, model.hy, 0.1, 0.04, 0.5, 0.5) + mu;
    // for (int j = 150; j < 175; j++)
    // {
    //     for (int i = 50; i < 150; i++)
    //     {
    //         eps(i,j) = 0.5;
    //     }
    // }
    
    MatrixXd vx = MatrixXd::Zero(model.Nx, model.Ny);
    MatrixXd vy = MatrixXd::Zero(model.Nx, model.Ny);
    MatrixXd pp = MatrixXd::Zero(model.Nx, model.Ny);
    model.SetParameters(eps, mu);
    model.SetInit(vx, vy, pp);

    // source
    // vector<Coor2D> c;
    // c.reserve(201);
    // for (int i = 51; i < 151; i++)
    // {
    //     Coor2D ccc(0,i);
    //     c.emplace_back(ccc);
    // }
    Coor2D c1(50,100);
    Coor2D c2(150,100);
    Coor2D c3(151,51);
    Coor2D c4(151,151);
    VectorXd f1 = Sine1D(model.Nt, model.tau, 5);
    VectorXd f2 = Sine1D(model.Nt, model.tau, 10);
    VectorXd f3 = Sine1D(model.Nt, model.tau, 15);
    VectorXd f4 = Sine1D(model.Nt, model.tau, 20);
    // vector<Coor2D> c{c1,c2,c3,c4};
    // vector<VectorXd> f{f1,f2,f3,f4};
    vector<Coor2D> c{c1,c2};
    vector<VectorXd> f{f2,f3};
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
    model.SaveAllRes(".\\data\\solData.txt", 10);

}