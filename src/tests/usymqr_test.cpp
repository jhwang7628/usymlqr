#include <memory>
#include <fstream>
#include <iostream>
#include "macros.h"
#include "usym_linear_solver.hpp"
#include "sparse_matrix.hpp"
#include "usym_tridiag.hpp"
#include "Eigen/Dense"
//##############################################################################
static const int N = 20;
using T = double; 
using USYM_Vector = typename Sparse_Matrix<T,N>::USYM_Vector;
using USYM_VectorX= typename Sparse_Matrix<T,N>::USYM_VectorX;
using USYM_Matrix = typename Sparse_Matrix<T,N>::USYM_Matrix; 
using USYM_MatrixX= typename Sparse_Matrix<T,N>::USYM_MatrixX; 
//##############################################################################
void CheckError(const USYM_Matrix &matA,  
                const USYM_Matrix &matP,  
                const USYM_Matrix &matQ, 
                const USYM_Matrix &matT,
                const int j) 
{
    using Matrix = typename Sparse_Matrix<T,N>::USYM_MatrixX; 
    using Vector = typename Sparse_Matrix<T,N>::USYM_Vector; 

    const Matrix Q_j = matQ.block(0,0,matQ.rows(),j); 
    const Matrix P_j = matP.block(0,0,matP.rows(),j); 
    const Eigen::MatrixXd T_j = matT.block(0,0,j,j); 
    const T beta_jp1  = matT(j  , j-1); 
    const T gamma_jp1 = matT(j-1, j  ); 
    const Vector p_jp1 = matP.col(j);
    const Vector q_jp1 = matQ.col(j);

    // check AQ_j = P_j S_j = P_j T_j +  beta_{j+1} p_{j+1} e_j^T
    Matrix AQ = matA * Q_j; 
    Matrix PT = P_j  * T_j; 
    Matrix PS = PT; PS.col(j-1) += beta_jp1*p_jp1;
    Matrix diff_1 = AQ - PS; 
    T res_1 = diff_1.norm(); 

    // check A^T P_j = Q_j S_j = Q_j T_j^T +  gamma_{j+1} q_{j+1} e_j^T
    Matrix AP = matA.transpose() * P_j; 
    Matrix QT = Q_j  * T_j.transpose(); 
    Matrix QS = QT; QS.col(j-1) += gamma_jp1*q_jp1;
    Matrix diff_2 = AP - QS; 
    T res_2 = diff_2.norm(); 

    std::cout << "***********************************\n";
    std::cout << "CHECK TRIDIAG ERROR: \n"; 
    std::cout << " Error for EQ2.6a = " << res_1 << "\n";
    std::cout << " Error for EQ2.6b = " << res_2 << "\n";
    std::cout << "***********************************\n";
}

//##############################################################################
void Basics() 
{
    std::cout << "\n========== Basics ===========\n"; 
    Eigen::Matrix3d m;
    m << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    Sparse_Matrix<double,10> sm; 
}
//##############################################################################
void Tridiagonalization(const USYM_Vector &b, const USYM_Vector &c,
                        const USYM_Matrix &matA, 
                              USYM_Matrix &matP, 
                              USYM_Matrix &matQ, 
                              USYM_Matrix &matT)
{
    auto sm = std::make_shared<Sparse_Matrix<T,N>>(matA); 
    USYM_Tridiag<T,N> tridiag(sm); 

    // test of complete tridiagonalization
    matP.setZero();
    matQ.setZero();
    matT.setZero(); 

    USYM_Vector p1,q1,p2,q2,p3,q3; 
    tridiag.Set_b(b); 
    tridiag.Set_c(c); 
    tridiag.InitialStep(p1, q1,
                        p2, q2, 
                        matT(0,0), matT(1,0), matT(0,1)); 
    matP.col(0) = p1; 
    matP.col(1) = p2; 
    matQ.col(0) = q1; 
    matQ.col(1) = q2; 
    for (int ii=1; ii<N; ++ii) // starting from i=2 loop
    {
        p1 = matP.col(ii-1); 
        p2 = matP.col(ii  ); 
        q1 = matQ.col(ii-1); 
        q2 = matQ.col(ii  ); 
        T tmp[2]; 
        bool finish = tridiag.Step(p1, q1,
                                   p2, q2,
                                   matT(ii,ii-1), matT(ii-1,ii),
                                   p3, q3,
                                   matT(ii,ii), tmp[0], tmp[1]
                                   ); 
        if (finish || ii==N-1)
            break; 

        matP.col(ii+1) = p3; 
        matQ.col(ii+1) = q3; 
        matT(ii+1,ii)= tmp[0];
        matT(ii,ii+1)= tmp[1]; 
    }
}
//##############################################################################
void Tridiagonalization()
{
    std::cout << "\n========== Tridiagonalization ===========\n"; 
    USYM_Matrix matA, matP, matQ, matT; 
    matA = USYM_Matrix::Random();
    USYM_Vector b = USYM_Vector::Random(); 
    USYM_Vector c = USYM_Vector::Random(); 
    Tridiagonalization(b, c, matA, matP, matQ, matT);
    CheckError(matA, matP, matQ, matT, N-1);
}
//##############################################################################
void Naive_Linear_Solve()
{
    std::cout << "\n========== Naive Linear Solve ===========\n"; 
    // declare
    USYM_Matrix matA, matP, matQ, matT; 
    USYM_Vector b, x0, r; 

    // initialize and compute tridiagonalization
    matA = USYM_Matrix::Random();
    b = USYM_Vector::Random();
    x0.setZero(); 
    r = (b - matA*x0);
    USYM_Vector init_b = r; 
    USYM_Vector init_c = USYM_Vector::Random(); 
    Tridiagonalization(init_b, init_c, matA, matP, matQ, matT);

    // naive solve with LQ
    const T beta_1 = r.norm(); 
    int jj = N;
    //for (int jj=2; jj<=N; ++jj)
    {
        Eigen::MatrixXd T_j = matT.block(0,0,jj,jj); 
        Eigen::VectorXd b_j = Eigen::VectorXd::Zero(jj); b_j[0] = beta_1; 
        // use full-pivot LU to solve T_j h_j_cg = beta_1 e_1
        Eigen::VectorXd h_j = T_j.fullPivLu().solve(b_j); 
        if(!(T_j*h_j).isApprox(b_j))
            std::cerr << "**ERROR** There is no solution.\n";

        if (jj<N)
        {
            T r_estimated = matT(jj,jj-1)*std::abs(h_j(jj-1)); 
            std::cout << "r_est = " << r_estimated << std::endl;
        }

        Eigen::VectorXd x_j_cg = x0 + matQ.block(0,0,matQ.rows(),jj)*h_j;
        r = (b - matA*x_j_cg); 
        std::cout << "\nstep " << jj << " has residual = " << r.norm() << std::endl; 
        std::cout << "T_j h_j - beta_1 e_1 = " << (T_j*h_j - b_j).norm() << std::endl;
        std::cout << "x_j_cg=" << x_j_cg.transpose() << std::endl; 
    }
    std::cout << "x_*   =" << matA.fullPivLu().solve(b).transpose() << std::endl;
}
//##############################################################################
void Linear_Solve()
{
    std::cout << "\n========== Linear Solve ===========\n"; 
    const int N = 5; 
    using T = double; 
    using USYM_Vector = typename Sparse_Matrix<double,N>::USYM_Vector;
    using USYM_VectorX= typename Sparse_Matrix<double,N>::USYM_VectorX;
    using USYM_Matrix = typename Sparse_Matrix<double,N>::USYM_Matrix; 
    using USYM_MatrixX= typename Sparse_Matrix<double,N>::USYM_MatrixX; 
    // initialize A and b
    USYM_Matrix A;
    A << 1, 0, 0, 1, 0,
         0, 1, 0, 1, 0,
         0, 0, 1, 0, 0,
         0, 0, 0, 1, 1,
         0, 0, 0, 0, 1; 
    USYM_Vector b; 
    b << 1, 2, 3, 4, 5; 
    // initialize solver with x0
    auto sparse_matrix = std::make_shared<Sparse_Matrix<T,N>>(A); 
    USYM_Linear_Solver<double,N> solver(sparse_matrix,b); 
    USYM_Vector x0 = USYM_Vector::Zero(); 
    solver.Initialize(x0); 
    for (int ii=0; ii<N; ++ii)
        solver.Step();
}
//##############################################################################
int main() {
    Basics();
    Tridiagonalization();
    Naive_Linear_Solve();
    Linear_Solve();
}
//##############################################################################
