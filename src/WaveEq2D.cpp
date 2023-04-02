#include "WaveEq2D.h"

// model class
ModelWaveEq2D::ModelWaveEq2D() {
    Nx = 101;
    Ny = 101;
    Nt = 201;
    hx = 0.01;
    hy = 0.01;
    tau = 0.001;
    C = MatrixXd::Ones(Nx, Ny);
    C2 = C.cwiseProduct(C);
    U_init = MatrixXd::Zero(Nx, Ny);
}

// member functions
void ModelWaveEq2D::SetSpatial(uint16_t set_Nx, uint16_t set_Ny, const double& set_hx, const double& set_hy) {
    Nx = set_Nx;
    Ny = set_Ny;
    hx = set_hx;
    hy = set_hy;
}

void ModelWaveEq2D::SetTime(uint16_t set_Nt, const double& set_tau) {
    Nt = set_Nt;
    tau = set_tau;
}

void ModelWaveEq2D::SetC(const MatrixXd& set_C) {
    C = set_C;
    C2 = C.cwiseProduct(C);
}

void ModelWaveEq2D::SetInit(const MatrixXd& set_U_init) {
    U_init = set_U_init;
}

void ModelWaveEq2D::SetSource(const PointSource2D& set_source) {
    source = set_source;
}

void ModelWaveEq2D::PrintInfo() const{
    std::cout << "Model size is " << Nx << " by " << Ny << std::endl;
    std::cout << "hx = " << hx << ", hy = " << hy << std::endl;
    std::cout << "Maximum velocity is: " << C.maxCoeff() << std::endl;
}

// wave equation class
WaveEq2D::WaveEq2D() {
    LU = MatrixXd::Zero(Nx, Ny);
    U0 = DiffGrid2D(Nx, Ny, hx, hy);
    U1 = DiffGrid2D(Nx, Ny, hx, hy);
    U2 = DiffGrid2D(Nx, Ny, hx, hy);
    source = PointSource2D();
    Record = false;
}

void WaveEq2D::SetRecord(bool r) {
    Record = r;
}

// spatial differential operator:
// LU = C^2 (U_xx + U_yy)
void WaveEq2D::EvalLU(DiffGrid2D& U) {
    LU = C2.cwiseProduct(U.dxx() + U.dyy());
}

void WaveEq2D::TimeUpdate(uint16_t idx_t) {
    if (source.GetSourceNum() == 0) {

        EvalLU(U1);
        U2.SetValue(2*U1.GetValue() - U0.GetValue() + tau*tau*LU);
        U0.SetValue(U1.GetValue());
        U1.SetValue(U2.GetValue());

    } else {

        EvalLU(U1);
        U2.SetValue(2*U1.GetValue() - U0.GetValue() + tau*tau*LU);
        for (int i = 0; i < source.GetSourceNum(); i++)
        {
            U2.AddOneEntry(source.GetCoor(i).x, source.GetCoor(i).y, tau*tau*source.GetValue(i, idx_t));
        }
        U0.SetValue(U1.GetValue());
        U1.SetValue(U2.GetValue());
    }
}

void WaveEq2D::Solve() {

    LU = MatrixXd::Zero(Nx, Ny);
    U0 = DiffGrid2D(U_init, hx, hy);
    U1 = DiffGrid2D(Nx, Ny, hx, hy);
    U2 = DiffGrid2D(Nx, Ny, hx, hy);

    if (Record == false) {
        for (int idx_t = 0; idx_t < Nt; idx_t++)
        {
            TimeUpdate(idx_t);
        }
    } else {
        SolData.reserve(Nt);
        for (int idx_t = 0; idx_t < Nt; idx_t++)
        {
            TimeUpdate(idx_t);
            SolData.emplace_back(U0.GetValue());
        }
    }
}

void WaveEq2D::SaveSol(std::string savename) {
    std::fstream file;
    file.open(savename, std::ios_base::out);
    for (auto e : U0.GetValue().reshaped() ) {
        file << e << " ";
    }
    file.close();
}

void WaveEq2D::SaveSolData(std::string savename, uint16_t step_t) {
    std::fstream file;
    file.open(savename, std::ios_base::out);
    for (int i = 0; i < Nt; i+=step_t)
    {
        for (auto e : SolData[i].reshaped()) {
            file << e << " ";
        }
        file << std::endl;
    }
    file.close();
}