#include "ConvDiffEq2D.h"

// constructor
ConvDiffEq2D::ConvDiffEq2D() {
    Record = false;

    D = MatrixXf::Ones(Nx, Ny);
    Vx = MatrixXf::Ones(Nx, Ny);
    Vy = MatrixXf::Ones(Nx, Ny);
    U_init = MatrixXf::Zero(Nx, Ny);

    source = PointSource2D();

    U0 = DiffGrid2D(Nx, Ny, hx, hy);
    U1 = DiffGrid2D(Nx, Ny, hx, hy);
    DdxU = DiffGrid2D(Nx, Ny, hx, hy);
    DdyU = DiffGrid2D(Nx, Ny, hx, hy);
    VxU = DiffGrid2D(Nx, Ny, hx, hy);
    VyU = DiffGrid2D(Nx, Ny, hx, hy);
}

// member functions
void ConvDiffEq2D::SetParameters(const MatrixXf& set_D, const MatrixXf set_Vx, const MatrixXf set_Vy){
    D = set_D;
    Vx = set_Vx;
    Vy = set_Vy;
}

void ConvDiffEq2D::SetInit(const MatrixXf& set_U_init) {
    U_init = set_U_init;
}

void ConvDiffEq2D::SetSource(const PointSource2D& set_source) {
    source = set_source;
}

void ConvDiffEq2D::SetRecord(bool r) {
    Record = r;
}

// solver

void ConvDiffEq2D::TimeUpdate(int idx_t) {

    if (source.GetSourceNum() == 0) {

        DdxU.SetValue(D.cwiseProduct(U0.dx_b()));
        DdyU.SetValue(D.cwiseProduct(U0.dy_b()));
        VxU.SetValue(Vx.cwiseProduct(U0.GetValue()));
        VyU.SetValue(Vy.cwiseProduct(U0.GetValue()));

        U1.SetValue(
            U0.GetValue()
            + tau * (DdxU.dx_f() + DdyU.dy_f() + D.cwiseProduct(U0.dxx() + U0.dyy()))
            - tau * (VxU.dx_f() + VyU.dy_f() + Vx.cwiseProduct(U0.dx_f()) + Vy.cwiseProduct(U0.dy_f()))
        );

        U0.SetValue(U1.GetValue());

    } else {

        DdxU.SetValue(D.cwiseProduct(U0.dx_f()));
        DdyU.SetValue(D.cwiseProduct(U0.dy_f()));
        VxU.SetValue(Vx.cwiseProduct(U0.GetValue()));
        VyU.SetValue(Vy.cwiseProduct(U0.GetValue()));

        U1.SetValue(
            U0.GetValue()
            + tau * (DdxU.dx_b() + DdyU.dy_b() + D.cwiseProduct(U0.dxx() + U0.dyy()))
            - tau * (VxU.dx_b() + VyU.dy_b() + Vx.cwiseProduct(U0.dx_f()) + Vy.cwiseProduct(U0.dy_f()))
        );

        // source
        for (int i = 0; i < source.GetSourceNum(); i++)
        {
            // U1.AddOneEntry(source.GetCoor(i).x, source.GetCoor(i).y, tau*source.GetValue(i, idx_t));
            U1.SetOneEntry(source.GetCoor(i).x, source.GetCoor(i).y, source.GetValue(i, idx_t));
        }

        U0.SetValue(U1.GetValue());

    }
}

void ConvDiffEq2D::Solve() {
    ConvDiffEq2D();
    U0 = DiffGrid2D(U_init, hx, hy);
    U1 = DiffGrid2D(Nx, Ny, hx, hy);
    DdxU = DiffGrid2D(Nx, Ny, hx, hy);
    DdyU = DiffGrid2D(Nx, Ny, hx, hy);
    VxU = DiffGrid2D(Nx, Ny, hx, hy);
    VyU = DiffGrid2D(Nx, Ny, hx, hy);

    if (Record == false)
    {
        for (int idx_t = 0; idx_t < Nt; idx_t++)
        {
            TimeUpdate(idx_t);
        }
        Res = ExtractFromPML(U0.GetValue());
    } else {
        AllRes.reserve(Nt);
        for (int idx_t = 0; idx_t < Nt; idx_t++)
        {
            TimeUpdate(idx_t);
            AllRes.emplace_back(ExtractFromPML(U0.GetValue()));
        }
        Res = ExtractFromPML(U0.GetValue());
    }
}