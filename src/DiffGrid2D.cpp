#include "DiffGrid2D.h"

using Eigen::MatrixXf;

DiffGrid2D::DiffGrid2D() {
    hx = hy = 0.;
    U = MatrixXf();
    Ux = MatrixXf();
    Uy = MatrixXf();
    Uxx = MatrixXf();
    Uyy = MatrixXf();
}

DiffGrid2D::DiffGrid2D(int Nx, int Ny, double hx, double hy) : hx(hx), hy(hy)
{
    U = MatrixXf::Zero(Nx, Ny);
    Ux = MatrixXf::Zero(Nx, Ny);
    Uy = MatrixXf::Zero(Nx, Ny);
    Uxx = MatrixXf::Zero(Nx, Ny);
    Uyy = MatrixXf::Zero(Nx, Ny);
}

DiffGrid2D::DiffGrid2D(const MatrixXf& set_U, double hx, double hy) : hx(hx), hy(hy)
{
    U = set_U;
    Ux = MatrixXf::Zero(U.rows(), U.cols());
    Uy = MatrixXf::Zero(U.rows(), U.cols());
    Uxx = MatrixXf::Zero(U.rows(), U.cols());
    Uyy = MatrixXf::Zero(U.rows(), U.cols());
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

void DiffGrid2D::SetValue(const MatrixXf& set_U) {
    U = set_U;
}

void DiffGrid2D::Seth(const double set_hx, const double set_hy) {
    hx = set_hx;
    hy = set_hy;
}

const MatrixXf& DiffGrid2D::GetValue() {
    return U;
}

// differential opeartor, forward and backward
const MatrixXf& DiffGrid2D::dx_f() {
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

const MatrixXf& DiffGrid2D::dx_b() {
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

const MatrixXf& DiffGrid2D::dy_f() {
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

const MatrixXf& DiffGrid2D::dy_b() {
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

const MatrixXf& DiffGrid2D::dxx() {
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

const MatrixXf& DiffGrid2D::dyy() {
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

void Model2D::SetSpatial(int set_Nx, int set_Ny, const double& set_hx, const double& set_hy) {
    Nx = set_Nx;
    Ny = set_Ny;
    hx = set_hx;
    hy = set_hy;
}

void Model2D::SetTime(int set_Nt, const double& set_tau) {
    Nt = set_Nt;
    tau = set_tau;
}

void Model2D::SetPML(int set_N_pml, const double& set_coef_pml) {
    N_pml = set_N_pml;
    coef_pml = set_coef_pml;
}

void Model2D::BuildSigmaPML() {
    sigmaX = MatrixXf::Zero(Nx+2*N_pml, Ny+2*N_pml);
    sigmaY = MatrixXf::Zero(Nx+2*N_pml, Ny+2*N_pml);
    for (int i = 0; i < N_pml; i++)
    {
        sigmaX.row(N_pml-1-i).setConstant(i*i*coef_pml);
        sigmaX.row(N_pml+Nx+i).setConstant(i*i*coef_pml);

        sigmaY.col(N_pml-1-i).setConstant(i*i*coef_pml);
        sigmaY.col(N_pml+Ny+i).setConstant(i*i*coef_pml);
    }
}

MatrixXf Model2D::ExpandToPML(const MatrixXf& M) {
    MatrixXf M_pml = MatrixXf::Zero(Nx+2*N_pml, Ny+2*N_pml);
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

MatrixXf Model2D::ExtractFromPML(const MatrixXf& M_pml) {
    MatrixXf M = MatrixXf::Zero(Nx, Ny);
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
    std::cout << "Total time: " << Nt * tau << " s" << std::endl;
    if (N_pml != 0) {
        std::cout << "N_pml = " << N_pml << ", coef_pml = " << coef_pml << std::endl;
    }
}

// save: do not optimize with O1, O2, or O3. will bypass the fstream file operation.
#pragma GCC push_options
#pragma GCC optimize("O0")

void Model2D::SaveRes(std::string savename) {
    std::fstream file;
    file.open(savename, std::ios_base::out);
    for (auto e : Res.reshaped() ) {
        file << e << " ";
    }
    file.close();
}

void Model2D::SaveAllRes(std::string savename, int step_t) {
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

#pragma GCC pop_options