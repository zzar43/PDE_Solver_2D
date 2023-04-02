#include "DiffGrid2D.h"

using Eigen::MatrixXd;

DiffGrid2D::DiffGrid2D() {
    hx = hy = 0.;
    U = MatrixXd();
    Ux = MatrixXd();
    Uy = MatrixXd();
    Uxx = MatrixXd();
    Uyy = MatrixXd();
}

DiffGrid2D::DiffGrid2D(uint16_t Nx, uint16_t Ny, double hx, double hy) : hx(hx), hy(hy)
{
    U = MatrixXd::Zero(Nx, Ny);
    Ux = MatrixXd::Zero(Nx, Ny);
    Uy = MatrixXd::Zero(Nx, Ny);
    Uxx = MatrixXd::Zero(Nx, Ny);
    Uyy = MatrixXd::Zero(Nx, Ny);
}

DiffGrid2D::DiffGrid2D(const MatrixXd& set_U, double hx, double hy) : hx(hx), hy(hy)
{
    U = set_U;
    Ux = MatrixXd::Zero(U.rows(), U.cols());
    Uy = MatrixXd::Zero(U.rows(), U.cols());
    Uxx = MatrixXd::Zero(U.rows(), U.cols());
    Uyy = MatrixXd::Zero(U.rows(), U.cols());
}

DiffGrid2D::DiffGrid2D(const DiffGrid2D& set_DiffGrid2D) {
    hx = set_DiffGrid2D.hx;
    hy = set_DiffGrid2D.hy;
    U = set_DiffGrid2D.U;
    Ux = set_DiffGrid2D.Ux;
    Uy = set_DiffGrid2D.Uy;
    Uxx = set_DiffGrid2D.Uxx;
    Uyy = set_DiffGrid2D.Uyy;
}

DiffGrid2D::~DiffGrid2D() {}

void DiffGrid2D::SetValue(const MatrixXd& set_U) {
    U = set_U;
}

void DiffGrid2D::Seth(const double set_hx, const double set_hy) {
    hx = set_hx;
    hy = set_hy;
}

const MatrixXd& DiffGrid2D::GetValue() {
    return U;
}

// differential opeartor, forward and backward
const MatrixXd& DiffGrid2D::dx_f() {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < U.cols()-1; j++)
    {
        for (int i = 1; i < U.rows()-1; i++)
        {
            Ux(i,j) = (U(i+1,j) - U(i,j)) / hx;
        }
    }
    return Ux;
}

const MatrixXd& DiffGrid2D::dx_b() {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < U.cols()-1; j++)
    {
        for (int i = 1; i < U.rows()-1; i++)
        {
            Ux(i,j) = (U(i,j) - U(i-1,j)) / hx;
        }
    }
    return Ux;
}

const MatrixXd& DiffGrid2D::dy_f() {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < U.cols()-1; j++)
    {
        for (int i = 1; i < U.rows()-1; i++)
        {
            Uy(i,j) = (U(i,j+1) - U(i,j)) / hy;
        }
    }
    return Uy;
}

const MatrixXd& DiffGrid2D::dy_b() {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < U.cols()-1; j++)
    {
        for (int i = 1; i < U.rows()-1; i++)
        {
            Uy(i,j) = (U(i,j) - U(i,j-1)) / hy;
        }
    }
    return Uy;
}

const MatrixXd& DiffGrid2D::dxx() {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < U.cols()-1; j++)
    {
        for (int i = 1; i < U.rows()-1; i++)
        {
            Uxx(i,j) = (U(i+1,j) - 2*U(i,j) + U(i-1,j)) / hx / hx;
        }
    }
    return Uxx;
}

const MatrixXd& DiffGrid2D::dyy() {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < U.cols()-1; j++)
    {
        for (int i = 1; i < U.rows()-1; i++)
        {
            Uyy(i,j) = (U(i,j+1) - 2*U(i,j) + U(i,j-1)) / hy / hy;
        }
    }
    return Uyy;
}

void DiffGrid2D::PrintShape() {
    std::cout << " Nx = " << U.rows() << " Ny = " << U.cols();
    std::cout << " hx = " << hx << " hy = " << hy << std::endl;
}

// a general 2D model class for defining basic parameters
Model2D::Model2D() {
    Nx = 101;
    Ny = 101;
    N_pml = 0;
    Nt = 201;
    hx = 0.01;
    hy = 0.01;
    tau = 0.001;
    coef_pml = 1.;
}

void Model2D::SetSpatial(uint16_t set_Nx, uint16_t set_Ny, const double& set_hx, const double& set_hy) {
    Nx = set_Nx;
    Ny = set_Ny;
    hx = set_hx;
    hy = set_hy;
}

void Model2D::SetTime(uint16_t set_Nt, const double& set_tau) {
    Nt = set_Nt;
    tau = set_tau;
}

void Model2D::SetPML(uint16_t set_N_pml, const double& set_coef_pml) {
    N_pml = set_N_pml;
    coef_pml = set_coef_pml;
}

void Model2D::BuildSigmaPML() {
    sigmaX = MatrixXd::Zero(Nx+2*N_pml, Ny+2*N_pml);
    sigmaY = MatrixXd::Zero(Nx+2*N_pml, Ny+2*N_pml);
    for (int i = 0; i < N_pml; i++)
    {
        sigmaX.row(N_pml-1-i).setConstant(i*i*coef_pml);
        sigmaX.row(N_pml+Nx+i).setConstant(i*i*coef_pml);

        sigmaY.col(N_pml-1-i).setConstant(i*i*coef_pml);
        sigmaY.col(N_pml+Ny+i).setConstant(i*i*coef_pml);
    }
}

MatrixXd Model2D::ExpandToPML(const MatrixXd& M) {
    MatrixXd M_pml = MatrixXd::Zero(Nx+2*N_pml, Ny+2*N_pml);
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            M_pml(i+N_pml, j+N_pml) = M(i,j);
        }
    }
    for (int i = 0; i < N_pml; i++)
    {
        M_pml.row(i) = M_pml.row(N_pml);
        M_pml.row(i+N_pml+Nx) = M_pml.row(N_pml+Nx-1);
    }
    for (int i = 0; i < N_pml; i++)
    {
        M_pml.col(i) = M_pml.col(N_pml);
        M_pml.col(i+N_pml+Ny) = M_pml.col(N_pml+Ny-1);
    }
    return M_pml;
}

MatrixXd Model2D::ExtractFromPML(const MatrixXd& M_pml) {
    MatrixXd M = MatrixXd::Zero(Nx, Ny);
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            M(i,j) = M_pml(i+N_pml, j+N_pml);
        }
    }
    return M;
}

// print
void Model2D::PrintInfo() {
    std::cout << "Model information: " << std::endl;
    std::cout << "Nx = " << Nx << ", Ny = " << Ny << std::endl;
    std::cout << "hx = " << hx << ", hy = " << hy << std::endl;
    std::cout << "Nt = " << Nt << ", tau = " << tau << std::endl;
    std::cout << "N_pml = " << N_pml << ", coef_pml = " << coef_pml << std::endl;
}

// save
void Model2D::SaveRes(std::string savename) {
    std::fstream file;
    file.open(savename, std::ios_base::out);
    for (auto e : Res.reshaped() ) {
        file << e << " ";
    }
    file.close();
}

void Model2D::SaveAllRes(std::string savename, uint16_t step_t) {
    std::fstream file;
    file.open(savename, std::ios_base::out);
    for (int i = 0; i < Nt; i+=step_t)
    {
        for (auto e : AllRes[i].reshaped()) {
            file << e << " ";
        }
        file << std::endl;
    }
    file.close();
}