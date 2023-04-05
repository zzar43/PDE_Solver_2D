#include "WaveEq2D.h"

// wave equation class
// constructor
WaveEq2D::WaveEq2D() {
    Record = false;

    LU = MatrixXf::Zero(Nx, Ny);
    C = MatrixXf::Ones(Nx, Ny);
    C2 = MatrixXf::Ones(Nx, Ny);
    U_init = MatrixXf::Zero(Nx, Ny);

    source = PointSource2D();

    U0 = DiffGrid2D(Nx, Ny, hx, hy);
    U1 = DiffGrid2D(Nx, Ny, hx, hy);
    U2 = DiffGrid2D(Nx, Ny, hx, hy);
}

void WaveEq2D::SetParameters(const MatrixXf& set_C) {
    C = set_C;
    C2 = C.cwiseProduct(C);
}

void WaveEq2D::SetInit(const MatrixXf& set_U_init) {
    U_init = set_U_init;
    U0.SetValue(U_init);
    U1.SetValue(U_init);
}

void WaveEq2D::SetSource(const PointSource2D& set_source) {
    source = set_source;
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
            // U2.SetOneEntry(source.GetCoor(i).x, source.GetCoor(i).y, tau*tau*source.GetValue(i, idx_t));
        }
        U0.SetValue(U1.GetValue());
        U1.SetValue(U2.GetValue());
    }
}

void WaveEq2D::Solve() {

    LU = MatrixXf::Zero(Nx, Ny);
    U0 = DiffGrid2D(U_init, hx, hy);
    U1 = DiffGrid2D(Nx, Ny, hx, hy);
    U2 = DiffGrid2D(Nx, Ny, hx, hy);

    if (Record == false) {
        for (int idx_t = 0; idx_t < Nt; idx_t++)
        {
            TimeUpdate(idx_t);
        }
        Res = U1.GetValue();
    } else {
        AllRes.reserve(Nt);
        for (int idx_t = 0; idx_t < Nt; idx_t++)
        {
            TimeUpdate(idx_t);
            AllRes.emplace_back(U1.GetValue());
        }
        Res = U1.GetValue();
    }
}
