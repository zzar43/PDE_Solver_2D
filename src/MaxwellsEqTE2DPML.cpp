#include "MaxwellsEqTE2DPML.h"


// constructor
MaxwellsEqTE2DPML::MaxwellsEqTE2DPML() {
    Record = false;

    Eps = MatrixXd::Ones(Nx+2*N_pml, Ny+2*N_pml);
    Mu = MatrixXd::Ones(Nx+2*N_pml, Ny+2*N_pml);
    A = Mu.cwiseInverse();
    B = Eps.cwiseInverse();

    Ex_init = MatrixXd::Zero(Nx+2*N_pml, Ny+2*N_pml);
    Ey_init = MatrixXd::Zero(Nx+2*N_pml, Ny+2*N_pml);
    Hz_init = MatrixXd::Zero(Nx+2*N_pml, Ny+2*N_pml);

    source = PointSource2D();

    Ex0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Ex1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Ey0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Ey1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Hz0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Hz1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    P0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    P1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Q0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Q1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    R0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    R1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
}

// member functions
void MaxwellsEqTE2DPML::SetParameters(const MatrixXd& set_Eps, const MatrixXd& set_Mu) {
    Eps = ExpandToPML(set_Eps);
    Mu = ExpandToPML(set_Mu);
    A = Mu.cwiseInverse();
    B = Eps.cwiseInverse();
}

void MaxwellsEqTE2DPML::SetInit(const MatrixXd& set_Vx_init, const MatrixXd& set_Vy_init, const MatrixXd& set_P_init) {
    Ex_init = ExpandToPML(set_Vx_init);
    Ey_init = ExpandToPML(set_Vy_init);
    Hz_init = ExpandToPML(set_P_init);
}

void MaxwellsEqTE2DPML::SetSource(const PointSource2D& set_source) {
    source = set_source;
    source.UpdatePML(N_pml);
}

void MaxwellsEqTE2DPML::SetRecord(bool r) {
    Record = r;
}

// solve
void MaxwellsEqTE2DPML::TimeUpdate(uint16_t idx_t) {

    if (source.GetSourceNum() == 0) {

        Ex1.SetValue(Ex0.GetValue() + tau*B.cwiseProduct(Hz0.dy_f()) - tau*sigmaY.cwiseProduct(Ex0.GetValue()));

        Ey1.SetValue(Ey0.GetValue() + -1*tau*B.cwiseProduct(Hz0.dx_f()) - tau*sigmaX.cwiseProduct(Ey0.GetValue()));

        Hz1.SetValue(Hz0.GetValue()
            - tau*Hz0.GetValue().cwiseProduct(sigmaX+sigmaY)
            + tau*sigmaX.cwiseProduct(sigmaY.cwiseProduct(P0.GetValue()))
            + tau*A.cwiseProduct(
                Ex0.dy_b() - Ey0.dx_b() - sigmaX.cwiseProduct(Q0.GetValue()) + sigmaY.cwiseProduct(R0.GetValue())
            ));

        P1.SetValue(P0.GetValue() - tau * Hz0.GetValue());
        Q1.SetValue(Q0.GetValue() - tau * Ex0.dy_f());
        R1.SetValue(R0.GetValue() - tau * Ey0.dx_f());

        Ex0.SetValue(Ex1.GetValue());
        Ey0.SetValue(Ey1.GetValue());
        Hz0.SetValue(Hz1.GetValue());
        P0.SetValue(P1.GetValue());
        Q0.SetValue(Q1.GetValue());
        R0.SetValue(R1.GetValue());

    } else {

        Ex1.SetValue(Ex0.GetValue() + tau*B.cwiseProduct(Hz0.dy_f()) - tau*sigmaY.cwiseProduct(Ex0.GetValue()));

        Ey1.SetValue(Ey0.GetValue() + -1*tau*B.cwiseProduct(Hz0.dx_f()) - tau*sigmaX.cwiseProduct(Ey0.GetValue()));

        Hz1.SetValue(Hz0.GetValue()
            - tau*Hz0.GetValue().cwiseProduct(sigmaX+sigmaY)
            + tau*sigmaX.cwiseProduct(sigmaY.cwiseProduct(P0.GetValue()))
            + tau*A.cwiseProduct(
                Ex1.dy_b() - Ey1.dx_b() + sigmaX.cwiseProduct(Q0.GetValue()) - sigmaY.cwiseProduct(R0.GetValue())
            ));

        // source
        for (int i = 0; i < source.GetSourceNum(); i++)
        {
            Hz1.AddOneEntry(source.GetCoor(i).x, source.GetCoor(i).y, tau*source.GetValue(i, idx_t));
            // Hz1.SetOneEntry(source.GetCoor(i).x, source.GetCoor(i).y, tau*source.GetValue(i, idx_t));
        }

        P1.SetValue(P0.GetValue() - tau * Hz0.GetValue());
        Q1.SetValue(Q0.GetValue() - tau * Ex0.dy_b());
        R1.SetValue(R0.GetValue() - tau * Ey0.dx_b());

        Ex0.SetValue(Ex1.GetValue());
        Ey0.SetValue(Ey1.GetValue());
        Hz0.SetValue(Hz1.GetValue());
        P0.SetValue(P1.GetValue());
        Q0.SetValue(Q1.GetValue());
        R0.SetValue(R1.GetValue());

    }
}

void MaxwellsEqTE2DPML::Solve() {
    MaxwellsEqTE2DPML();
    Ex0 = DiffGrid2D(Ex_init, hx, hy);
    Ex1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Ey0 = DiffGrid2D(Ey_init, hx, hy);
    Ey1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    Hz0 = DiffGrid2D(Hz_init, hx, hy);
    Hz1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    P1 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
    P0 = DiffGrid2D(Nx+2*N_pml, Ny+2*N_pml, hx, hy);
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
        Res = ExtractFromPML(Hz0.GetValue());
    } else {
        AllRes.reserve(Nt);
        for (int idx_t = 0; idx_t < Nt; idx_t++)
        {
            TimeUpdate(idx_t);
            AllRes.emplace_back(ExtractFromPML(Hz0.GetValue()));
        }
        Res = ExtractFromPML(Hz0.GetValue());
    }
}