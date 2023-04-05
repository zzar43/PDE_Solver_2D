#include <iostream>
#include <chrono>
#include <fstream>

#include "Eigen/Dense"
#include "DiffGrid2D.h"
// #include "WaveEq2D.h"
// #include "AcousticEq2DPML.h"
#include "MaxwellsEqTE2DPML.h"
#include "ConvDiffEq2D.h"

using Eigen::MatrixXf;
using Eigen::VectorXf;

using namespace std;
using namespace chrono; 

void NewSave(std::vector<MatrixXf> AllRes, std::string savename, int step_t) {
    std::fstream file;
    file.open(savename, std::ios_base::out);
    for (int i = 0; i < AllRes.size(); i+=step_t)
    {
        std::cout << "saving." << i << std::endl;
        for (auto e : AllRes[i].reshaped()) {
            file << e << " ";
        }
        file << std::endl;
    }
    file.close();
}

int main() {

    int sw = 1;

    // Convection-diffusion equation
    if (sw == 0) {
        ConvDiffEq2D model1;

        // spatial and time
        model1.SetSpatial(201,201,0.005,0.005);
        model1.SetTime(3001,0.0001);

        // model parameter
        MatrixXf vx = 0.5*MatrixXf::Ones(model1.Nx, model1.Ny);
        MatrixXf vy = 0*MatrixXf::Ones(model1.Nx, model1.Ny);
        MatrixXf D = 0.02*MatrixXf::Ones(model1.Nx, model1.Ny);
        for (int j = 120; j < 160; j++)
        {
            for (int i = 80; i < 120; i++)
            {
                D(j,i) = 0.01;
                vx(j,i) = 0.1;
            }
        }
        MatrixXf U_init = MatrixXf::Zero(model1.Nx, model1.Ny);
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
        VectorXf f1 = 0.01*VectorXf::Ones(model1.Nt);
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
    }

    // Maxwell's equations in 2D, transverse electric (TE) mode
    if (sw == 1) {

        MaxwellsEqTE2DPML model2;

        // set spatial, time and pml
        model2.SetSpatial(201,201,0.005,0.005);
        model2.SetTime(1001,0.001);
        model2.SetPML(30,0.1);

        // model parameter
        MatrixXf mu = MatrixXf::Ones(model2.Nx, model2.Ny);
        MatrixXf eps = 1*Gaussian2D(model2.Nx, model2.Ny, model2.hx, model2.hy, 0.1, 0.04, 0.5, 0.5) + mu;
        MatrixXf ex = MatrixXf::Zero(model2.Nx, model2.Ny);
        MatrixXf ey = MatrixXf::Zero(model2.Nx, model2.Ny);
        MatrixXf hz = MatrixXf::Zero(model2.Nx, model2.Ny);
        model2.SetParameters(eps,mu);
        model2.SetInit(ex,ey,hz);

        // source
        Coor2D c21(50,100);
        Coor2D c22(150,100);
        VectorXf f21 = Sine1D(model2.Nt, model2.tau, 10);
        VectorXf f22 = Sine1D(model2.Nt, model2.tau, 15);
        vector<Coor2D> c2{c21,c22};
        vector<VectorXf> f2{f21,f22};
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
    }


}