#include <iostream>
#include <chrono>

#include "Eigen/Dense"
#include "DiffGrid2D.h"
#include "WaveEq2D.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using namespace chrono; 

int main() {

    WaveEq2D model;

    model.SetSpatial(201,201,0.005,0.005);
    model.SetTime(2001,0.0005);
    
    MatrixXd cc = MatrixXd::Ones(model.Nx, model.Ny);
    MatrixXd uu = MatrixXd::Zero(model.Nx, model.Ny);
    model.SetC(cc);
    model.SetInit(uu);

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
    PointSource2D pp(c, f);
    model.SetSource(pp);

    // solve
    model.SetRecord(true);

    auto time1 = high_resolution_clock::now();
    model.Solve();
    auto time2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(time2 - time1);
    cout << "Computation time: " << duration.count() << " milliseconds"<< endl;

    // record
    model.SaveSol(".\\data\\sol.txt");
    model.SaveSolData(".\\data\\solData.txt", 100);

}