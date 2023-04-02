#include "AcousticEq2DPML.h"

// acoustic wave equation with PML class

// constructor
AcousticEq2DPML::AcousticEq2DPML() {
    Record = false;

    C = MatrixXd::Ones(Nx+2*N_pml, Ny+2*N_pml);
    Rho = MatrixXd::Ones(Nx+2*N_pml, Ny+2*N_pml);
    A = -1. * Rho.cwiseInverse();
    B = -1. * Rho.cwiseProduct(C.cwiseProduct(C));

    Vx_init = MatrixXd::Zero(Nx+2*N_pml, Ny+2*N_pml);
    Vy_init = MatrixXd::Zero(Nx+2*N_pml, Ny+2*N_pml);
    P_init = MatrixXd::Zero(Nx+2*N_pml, Ny+2*N_pml);

    source = PointSource2D();

    Vx0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Vx1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Vy0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Vy1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    P0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    P1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Q0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Q1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    R0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    R1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
}

// member functions
void AcousticEq2DPML::SetParameters(const MatrixXd& set_C, const MatrixXd& set_Rho) {
    C = ExpandToPML(set_C);
    Rho = ExpandToPML(set_Rho);
    A = -1. * Rho.cwiseInverse();
    B = -1. * Rho.cwiseProduct(C.cwiseProduct(C));
}

void AcousticEq2DPML::SetInit(const MatrixXd& set_Vx_init, const MatrixXd& set_Vy_init, const MatrixXd& set_P_init) {
    Vx_init = ExpandToPML(set_Vx_init);
    Vy_init = ExpandToPML(set_Vy_init);
    P_init = ExpandToPML(set_P_init);
}

void AcousticEq2DPML::SetSource(const PointSource2D& set_source) {
    source = set_source;
    source.UpdatePML(N_pml);
}

void AcousticEq2DPML::SetRecord(bool r) {
    Record = r;
}

// solve
void AcousticEq2DPML::TimeUpdate(uint16_t idx_t) {

    if (source.GetSourceNum() == 0) {

        Vx1.SetValue(Vx0.GetValue() + tau*A.cwiseProduct(P0.dx_f()) - tau*C.cwiseProduct(sigmaX.cwiseProduct(Vx0.GetValue())));
        Vy1.SetValue(Vy0.GetValue() + tau*A.cwiseProduct(P0.dy_f()) - tau*C.cwiseProduct(sigmaY.cwiseProduct(Vy0.GetValue())));
        Q1.SetValue(Q0.GetValue() + tau * C.cwiseProduct(Vy0.GetValue()));
        R1.SetValue(R0.GetValue() + tau * C.cwiseProduct(Vx0.GetValue()));
        P1.SetValue(P0.GetValue() 
            + tau * B.cwiseProduct(Vx1.dx_b() + Vy1.dy_b() + sigmaX.cwiseProduct(Q0.dy_b()) + sigmaY.cwiseProduct(R0.dx_b())) 
            - tau * C.cwiseProduct(P0.GetValue().cwiseProduct(sigmaX+sigmaY))
            );
        Vx0.SetValue(Vx1.GetValue());
        Vy0.SetValue(Vy1.GetValue());
        Q0.SetValue(Q1.GetValue());
        R0.SetValue(R1.GetValue());
        P0.SetValue(P1.GetValue());

    } else {

        Vx1.SetValue(Vx0.GetValue() + tau*A.cwiseProduct(P0.dx_f()) - tau*C.cwiseProduct(sigmaX.cwiseProduct(Vx0.GetValue())));
        Vy1.SetValue(Vy0.GetValue() + tau*A.cwiseProduct(P0.dy_f()) - tau*C.cwiseProduct(sigmaY.cwiseProduct(Vy0.GetValue())));
        Q1.SetValue(Q0.GetValue() + tau * C.cwiseProduct(Vy0.GetValue()));
        R1.SetValue(R0.GetValue() + tau * C.cwiseProduct(Vx0.GetValue()));
        P1.SetValue(P0.GetValue() 
            + tau * B.cwiseProduct(Vx1.dx_b() + Vy1.dy_b() + sigmaX.cwiseProduct(Q0.dy_b()) + sigmaY.cwiseProduct(R0.dx_b())) 
            - tau * C.cwiseProduct(P0.GetValue().cwiseProduct(sigmaX+sigmaY))
            );
        // source
        for (int i = 0; i < source.GetSourceNum(); i++)
        {
            P1.AddOneEntry(source.GetCoor(i).x, source.GetCoor(i).y, tau*source.GetValue(i, idx_t));
        }
        Vx0.SetValue(Vx1.GetValue());
        Vy0.SetValue(Vy1.GetValue());
        Q0.SetValue(Q1.GetValue());
        R0.SetValue(R1.GetValue());
        P0.SetValue(P1.GetValue());

    }
}

void AcousticEq2DPML::Solve() {
    AcousticEq2DPML();
    Vx0 = DiffGrid2D(Vx_init, hx, hy);
    Vx1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Vy0 = DiffGrid2D(Vy_init, hx, hy);
    Vy1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    P0 = DiffGrid2D(P_init, hx, hy);
    P1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Q1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Q0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    R1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    R0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    BuildSigmaPML();

    if (Record == false)
    {
        for (int idx_t = 0; idx_t < Nt; idx_t++)
        {
            TimeUpdate(idx_t);
        }
        Res = ExtractFromPML(P0.GetValue());
    } else {
        AllRes.reserve(Nt);
        for (int idx_t = 0; idx_t < Nt; idx_t++)
        {
            TimeUpdate(idx_t);
            AllRes.emplace_back(ExtractFromPML(P0.GetValue()));
        }
        Res = ExtractFromPML(P0.GetValue());
    }
}