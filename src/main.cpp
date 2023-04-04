#include <iostream>
#include <chrono>

#include "Eigen/Dense"
#include "DiffGrid2D.h"
// #include "WaveEq2D.h"
// #include "AcousticEq2DPML.h"
#include "MaxwellsEqTE2DPML.h"
#include "ConvDiffEq2D.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using namespace chrono; 


int main() {

    // ex1: convection-diffusion equation

    ConvDiffEq2D model1;

    // spatial and time
    model1.SetSpatial(201,201,0.005,0.005);
    model1.SetTime(3001,0.0001);

    // model parameter
    MatrixXd vx = 0.5*MatrixXd::Ones(model1.Nx, model1.Ny);
    MatrixXd vy = 0*MatrixXd::Ones(model1.Nx, model1.Ny);
    MatrixXd D = 0.02*MatrixXd::Ones(model1.Nx, model1.Ny);
    for (int j = 120; j < 160; j++)
    {
        for (int i = 80; i < 120; i++)
        {
            D(j,i) = 0.01;
            vx(j,i) = 0.1;
        }
    }
    MatrixXd U_init = MatrixXd::Zero(model1.Nx, model1.Ny);
    model1.SetParameters(D, vx, vy);
    model1.SetInit(U_init);

    // source
    vector<Coor2D> c1;
    c1.reserve(201);
    for (int i = 50; i < 150; i++)
    {
        Coor2D c_temp1(100,i);
        c1.emplace_back(c_temp1);
    }
    VectorXd f1 = 0.01*VectorXd::Ones(model1.Nt);
    PointSource2D ptsource1(c1, f1);
    model1.SetSource(ptsource1);

    // solve
    model1.SetRecord(true);
    cout << "Solving convection-diffusion equation." << endl;
    model1.PrintInfo();
    auto time1 = high_resolution_clock::now();
    model1.Solve();
    auto time2 = high_resolution_clock::now();
    auto duration1 = duration_cast<milliseconds>(time2 - time1);
    cout << "Computation time: " << duration1.count() << " milliseconds"<< endl;
    cout << "Done.\n" << endl;

    // save
    model1.SaveAllRes(".\\data\\solData_heat.txt", 50);

    // Maxwell's equations in 2D, transverse electric (TE) mode
    MaxwellsEqTE2DPML model2;

    // set spatial, time and pml
    model2.SetSpatial(201,201,0.005,0.005);
    model2.SetTime(1001,0.001);
    model2.SetPML(30,0.1);

    // model parameter
    MatrixXd mu = MatrixXd::Ones(model2.Nx, model2.Ny);
    MatrixXd eps = 1*Gaussian2D(model2.Nx, model2.Ny, model2.hx, model2.hy, 0.1, 0.04, 0.5, 0.5) + mu;
    MatrixXd ex = MatrixXd::Zero(model2.Nx, model2.Ny);
    MatrixXd ey = MatrixXd::Zero(model2.Nx, model2.Ny);
    MatrixXd hz = MatrixXd::Zero(model2.Nx, model2.Ny);
    model2.SetParameters(eps,mu);
    model2.SetInit(ex,ey,hz);

    // source
    Coor2D c21(50,100);
    Coor2D c22(150,100);
    VectorXd f21 = Sine1D(model2.Nt, model2.tau, 10);
    VectorXd f22 = Sine1D(model2.Nt, model2.tau, 15);
    vector<Coor2D> c2{c21,c22};
    vector<VectorXd> f2{f21,f22};
    PointSource2D ptsource2(c2, f2);
    model2.SetSource(ptsource2);

    // solve
    model2.SetRecord(true);
    cout << "Solving Maxwell's equations in 2D." << endl;
    model2.PrintInfo();
    auto time3 = high_resolution_clock::now();
    model2.Solve();
    auto time4 = high_resolution_clock::now();
    auto duration2 = duration_cast<milliseconds>(time4 - time3);
    cout << "Computation time: " << duration2.count() << " milliseconds"<< endl;
    cout << "Done.\n" << endl;

    // save
    model2.SaveAllRes(".\\data\\solData_wave.txt", 20);

    // // other

    // model.SetSpatial(201,201,0.005,0.005);
    // model.SetTime(3001,0.0001);
    // // model.SetPML(30,0.1);
    
    // // MatrixXd mu = MatrixXd::Ones(model.Nx, model.Ny);
    // // MatrixXd eps = MatrixXd::Ones(model.Nx, model.Ny);
    // // MatrixXd eps = 1*Gaussian2D(model.Nx, model.Ny, model.hx, model.hy, 0.1, 0.04, 0.5, 0.5) + mu;
    // // for (int j = 150; j < 175; j++)
    // // {
    // //     for (int i = 50; i < 150; i++)
    // //     {
    // //         eps(i,j) = 0.5;
    // //     }
    // // }
    
    // MatrixXd vx = 0.5*MatrixXd::Ones(model.Nx, model.Ny);
    // MatrixXd vy = 0*MatrixXd::Ones(model.Nx, model.Ny);
    // MatrixXd D = 0.02*MatrixXd::Ones(model.Nx, model.Ny);
    // for (int j = 120; j < 160; j++)
    // {
    //     for (int i = 80; i < 120; i++)
    //     {
    //         D(j,i) = 0.01;
    //         vx(j,i) = 0.1;
    //     }
    // }
    // MatrixXd U_init = MatrixXd::Zero(model.Nx, model.Ny);
    // MatrixXd U_init1 = 1*Gaussian2D(model.Nx, model.Ny, model.hx, model.hy, 0.1, 0.05, 0.3, 0.3);
    // MatrixXd U_init2 = 1*Gaussian2D(model.Nx, model.Ny, model.hx, model.hy, 0.05, 0.1, 0.7, 0.7);

    // model.SetParameters(D, vx, vy);
    // // model.SetInit(U_init1 + U_init2);
    // model.SetInit(U_init);

    // // source
    // vector<Coor2D> c;
    // c.reserve(201);
    // for (int i = 50; i < 150; i++)
    // {
    //     Coor2D ccc(100,i);
    //     c.emplace_back(ccc);
    // }
    // // VectorXd f = Gaussian1D(model.Nt, model.tau, 0.1, 0.1);
    // VectorXd f = 0.01*VectorXd::Ones(model.Nt);
    // // Coor2D c1(50,100);
    // // Coor2D c2(150,100);
    // // Coor2D c3(151,51);
    // // Coor2D c4(151,151);
    // // VectorXd f1 = Sine1D(model.Nt, model.tau, 5);
    // // VectorXd f2 = Sine1D(model.Nt, model.tau, 10);
    // // VectorXd f3 = Sine1D(model.Nt, model.tau, 15);
    // // VectorXd f4 = Sine1D(model.Nt, model.tau, 20);
    // // vector<Coor2D> c{c1,c2,c3,c4};
    // // vector<VectorXd> f{f1,f2,f3,f4};
    // // vector<Coor2D> c{c1,c2};
    // // vector<VectorXd> f{f2,f3};
    // PointSource2D ptsource(c, f);
    // model.SetSource(ptsource);



    // // solve
    // model.SetRecord(true);

    // auto time1 = high_resolution_clock::now();
    // model.Solve();
    // auto time2 = high_resolution_clock::now();
    // auto duration = duration_cast<milliseconds>(time2 - time1);
    // cout << "Computation time: " << duration.count() << " milliseconds"<< endl;

    // // record
    // model.SaveRes(".\\data\\sol.txt");
    // model.SaveAllRes(".\\data\\solData.txt", 50);

}