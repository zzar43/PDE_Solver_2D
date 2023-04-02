#include "AcousticEq2D.h"

// acoustic wave equation class
// constructor
AcousticEq2D::AcousticEq2D() {
    Record = false;

    C = MatrixXd::Ones(Nx, Ny);
    Rho = MatrixXd::Ones(Nx, Ny);
    A = -1. * Rho.cwiseInverse();
    B = -1. * Rho.cwiseProduct(C.cwiseProduct(C));

    Vx_init = MatrixXd::Zero(Nx, Ny);
    Vy_init = MatrixXd::Zero(Nx, Ny);
    P_init = MatrixXd::Zero(Nx, Ny);

    source = PointSource2D();

    Vx0 = DiffGrid2D(Nx, Ny, hx, hy);
    Vx1 = DiffGrid2D(Nx, Ny, hx, hy);
    Vy0 = DiffGrid2D(Nx, Ny, hx, hy);
    Vy1 = DiffGrid2D(Nx, Ny, hx, hy);
    P0 = DiffGrid2D(Nx, Ny, hx, hy);
    P1 = DiffGrid2D(Nx, Ny, hx, hy);
}

void AcousticEq2D::SetParameters(const MatrixXd& set_C, const MatrixXd& set_Rho) {
    C = set_C;
    Rho = set_Rho;
    A = -1. * Rho.cwiseInverse();
    B = -1. * Rho.cwiseProduct(C.cwiseProduct(C));
}

void AcousticEq2D::SetInit(const MatrixXd& set_Vx_init, const MatrixXd& set_Vy_init, const MatrixXd& set_P_init) {
    Vx_init = set_Vx_init;
    Vy_init = set_Vy_init;
    P_init = set_P_init;
}

void AcousticEq2D::SetSource(const PointSource2D& set_source) {
    source = set_source;
}

void AcousticEq2D::SetRecord(bool r) {
    Record = r;
}

void AcousticEq2D::TimeUpdate(uint16_t idx_t) {

    if (source.GetSourceNum() == 0) {

        Vx1.SetValue(Vx0.GetValue() + tau*A.cwiseProduct(P0.dx_f()));
        Vy1.SetValue(Vy0.GetValue() + tau*A.cwiseProduct(P0.dy_f()));
        P1.SetValue(P0.GetValue() + tau * B.cwiseProduct( Vx1.dx_b() + Vy1.dy_b() ));
        Vx0.SetValue(Vx1.GetValue());
        Vy0.SetValue(Vy1.GetValue());
        P0.SetValue(P1.GetValue());

    } else {

        Vx1.SetValue(Vx0.GetValue() + tau*A.cwiseProduct(P0.dx_f()));
        Vy1.SetValue(Vy0.GetValue() + tau*A.cwiseProduct(P0.dy_f()));
        P1.SetValue(P0.GetValue() + tau * B.cwiseProduct( Vx1.dx_b() + Vy1.dy_b() ));
        // source
        for (int i = 0; i < source.GetSourceNum(); i++)
        {
            P1.AddOneEntry(source.GetCoor(i).x, source.GetCoor(i).y, tau*source.GetValue(i, idx_t));
        }
        Vx0.SetValue(Vx1.GetValue());
        Vy0.SetValue(Vy1.GetValue());
        P0.SetValue(P1.GetValue());

    }
}

void AcousticEq2D::Solve() {
    AcousticEq2D();
    Vx0 = DiffGrid2D(Vx_init, hx, hy);
    Vx1 = DiffGrid2D(Nx, Ny, hx, hy);
    Vy0 = DiffGrid2D(Vy_init, hx, hy);
    Vy1 = DiffGrid2D(Nx, Ny, hx, hy);
    P0 = DiffGrid2D(P_init, hx, hy);
    P1 = DiffGrid2D(Nx, Ny, hx, hy);

    if (Record == false)
    {
        for (int idx_t = 0; idx_t < Nt; idx_t++)
        {
            TimeUpdate(idx_t);
        }
        Res = P0.GetValue();
    } else {
        AllRes.reserve(Nt);
        for (int idx_t = 0; idx_t < Nt; idx_t++)
        {
            TimeUpdate(idx_t);
            AllRes.emplace_back(P0.GetValue());
        }
        Res = P0.GetValue();
    }
}

// void AcousticEq2D::SaveSol(std::string savename) {
//     std::fstream file;
//     file.open(savename, std::ios_base::out);
//     for (auto e : P0.GetValue().reshaped() ) {
//         file << e << " ";
//     }
//     file.close();
// }

// void AcousticEq2D::SaveSolData(std::string savename, uint16_t step_t) {
//     std::fstream file;
//     file.open(savename, std::ios_base::out);
//     for (int i = 0; i < Nt; i+=step_t)
//     {
//         for (auto e : SolData[i].reshaped()) {
//             file << e << " ";
//         }
//         file << std::endl;
//     }
//     file.close();
// }
