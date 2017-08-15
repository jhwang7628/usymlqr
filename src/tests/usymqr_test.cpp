#include <memory>
#include <fstream>
#include <iostream>
#include "macros.h"
#include "usym_linear_solver.hpp"
#include "sparse_matrix.hpp"
#include "usym_tridiag.hpp"
#include "Eigen/Dense"
//##############################################################################
using T = double; 
using Vector = Eigen::Matrix<T,Eigen::Dynamic,1>; 
using Matrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>; 
//##############################################################################
template<class T_Vector, class T_Matrix>
void CheckError(const T_Matrix &matA,  
                const T_Matrix &matP,  
                const T_Matrix &matQ, 
                const T_Matrix &matT,
                const int j) 
{
    const T_Matrix Q_j = matQ.block(0,0,matQ.rows(),j); 
    const T_Matrix P_j = matP.block(0,0,matP.rows(),j); 
    const T_Matrix T_j = matT.block(0,0,j,j); 
    const T beta_jp1  = matT(j  , j-1); 
    const T gamma_jp1 = matT(j-1, j  ); 
    const T_Vector p_jp1 = matP.col(j);
    const T_Vector q_jp1 = matQ.col(j);

    // check AQ_j = P_j S_j = P_j T_j +  beta_{j+1} p_{j+1} e_j^T
    T_Matrix AQ = matA * Q_j; 
    T_Matrix PT = P_j  * T_j; 
    T_Matrix PS = PT; PS.col(j-1) += beta_jp1*p_jp1;
    T_Matrix diff_1 = AQ - PS; 
    T res_1 = diff_1.norm(); 

    // check A^T P_j = Q_j S_j = Q_j T_j^T +  gamma_{j+1} q_{j+1} e_j^T
    T_Matrix AP = matA.transpose() * P_j; 
    T_Matrix QT = Q_j  * T_j.transpose(); 
    T_Matrix QS = QT; QS.col(j-1) += gamma_jp1*q_jp1;
    T_Matrix diff_2 = AP - QS; 
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
}
//##############################################################################
template<class T_Vector, class T_Matrix>
void Tridiagonalization(const T_Vector &b, const T_Vector &c,
                        const std::shared_ptr<T_Matrix> &matA, 
                              T_Matrix &matP, 
                              T_Matrix &matQ, 
                              T_Matrix &matT)
{
    const int N = b.size();
    USYM_Tridiag<T,T_Vector,T_Matrix> tridiag(matA); 

    // test of complete tridiagonalization
    matP.setZero();
    matQ.setZero();
    matT.setZero(); 

    T_Vector p1,q1,p2,q2,p3,q3; 
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
template<class T_Vector, class T_Matrix>
void Tridiagonalization(const int N)
{
    std::cout << "\n========== Tridiagonalization ===========\n"; 
    std::shared_ptr<T_Matrix> matA(new T_Matrix(N,N));
    T_Matrix matP(N,N), matQ(N,N), matT(N,N); 
    (*matA) = T_Matrix::Random(N,N);
    T_Vector b = T_Vector::Random(N); 
    T_Vector c = T_Vector::Random(N); 
    Tridiagonalization<T_Vector,T_Matrix>(b, c, matA, matP, matQ, matT);
    CheckError<T_Vector,T_Matrix>(*matA, matP, matQ, matT, N-1);
}
//##############################################################################
template<class T_Vector, class T_Matrix>
T Naive_Linear_Solve(const int N)
{
    std::cout << "\n========== Naive Linear Solve ===========\n"; 
    // declare
    std::shared_ptr<T_Matrix> matA(new T_Matrix(N,N)); 
    T_Matrix matP(N,N), matQ(N,N), matT(N,N); 
    T_Vector b(N), x0(N), r(N); 

    // initialize and compute tridiagonalization
    (*matA) = T_Matrix::Random(N,N);
    b = T_Vector::Random(N);
    x0.setZero(); 
    r = (b - (*matA)*x0);
    T_Vector init_b = r; 
    T_Vector init_c = T_Vector::Random(N); 
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
        r = (b - (*matA)*x_j_cg); 
        std::cout << "\nstep " << jj << " has residual = " << r.norm() << std::endl; 
        std::cout << "T_j h_j - beta_1 e_1 = " << (T_j*h_j - b_j).norm() << std::endl;
        std::cout << "x_j_cg=" << x_j_cg.transpose() << std::endl; 
    }
    std::cout << "x_*   =" << matA->fullPivLu().solve(b).transpose() << std::endl;
    return r.norm();
}
//##############################################################################
template<class T_Vector, class T_Matrix>
void Linear_Solve(const int N)
{
    std::cout << "\n========== Linear Solve ===========\n"; 
    // initialize A and b
    std::shared_ptr<T_Matrix> A(new T_Matrix(N,N));
    (*A) << 1, 0, 0, 1, 0,
            0, 1, 0, 1, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 1,
            0, 0, 0, 0, 1; 
    T_Vector b(N); 
    b << 1, 2, 3, 4, 5; 
    // initialize solver with x0
    USYM_Linear_Solver<T,T_Vector,T_Matrix> solver(A,b); 
    T_Vector x0 = T_Vector::Zero(N); 
    solver.Initialize(x0); 
    for (int ii=0; ii<N; ++ii)
        solver.Step();
}
//##############################################################################
void Test_LossOrthogonality(const int maxN)
{
    std::cout << "\n========== TEST: Loss Orthogonality ===========\n"; 
    srand((unsigned int) time(0));
    for (int ii=0; ii<10; ++ii)
    {
        std::string filename("data/test_loss_ortho_" 
                             +std::to_string(ii)+".txt");
        std::ofstream ofs(filename);
        ofs << "N number_iterations residual\n";
        for (int N=2; N<maxN; ++N)
        {
            const T res = Naive_Linear_Solve<Vector,Matrix>(N); 
            ofs << N << " " << N << " " << res << std::endl;
        }
        ofs.close();
    }
}
//##############################################################################
int main() {
    const int N = 150;
    Basics();
    Tridiagonalization<Vector,Matrix>(N);
    Naive_Linear_Solve<Vector,Matrix>(N);
    Linear_Solve<Vector,Matrix>(5);

    // more complicated tests
    Test_LossOrthogonality(N);
}
//##############################################################################
