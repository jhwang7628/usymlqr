#include <memory>
#include <iostream>
#include "macros.h"
#include "sparse_matrix.hpp"
#include "usym_tridiag.hpp"
#include "Eigen/Dense"
//##############################################################################
void Basics() {
    Eigen::Matrix3d m;
    m << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    Sparse_Matrix<double,10> sm; 
}
//##############################################################################
void Tridiagonalization() {
    const int N = 5; 
    using USYM_Vector = typename Sparse_Matrix<double,N>::USYM_Vector;
    using USYM_VectorX= typename Sparse_Matrix<double,N>::USYM_VectorX;
    using USYM_Matrix = typename Sparse_Matrix<double,N>::USYM_Matrix; 
    using USYM_MatrixX= typename Sparse_Matrix<double,N>::USYM_MatrixX; 

    USYM_Matrix A = USYM_Matrix::Random();
    auto sm = std::make_shared<Sparse_Matrix<double,N>>(A); 
    USym_Tridiag<double,N> tridiag(sm); 

    // test of complete tridiagonalization
    USYM_Matrix P, Q, T; 
    P.setZero();
    Q.setZero();
    T.setZero(); 

    USYM_Vector p1,q1,p2,q2,p3,q3; 
    USYM_Vector b = USYM_Vector::Random(); 
    USYM_Vector c = USYM_Vector::Random(); 
    tridiag.InitialStep(b, c, 
                        p1, q1,
                        p2, q2, 
                        T(0,0), T(1,0), T(0,1)); 
    P.col(0) = p1; 
    P.col(1) = p2; 
    Q.col(0) = q1; 
    Q.col(1) = q2; 
    for (int ii=1; ii<N; ++ii)
    {
        p1 = P.col(ii-1); 
        p2 = P.col(ii  ); 
        q1 = Q.col(ii-1); 
        q2 = Q.col(ii  ); 
        double tmp[2]; 
        bool finish = tridiag.Step(p1, q1,
                                   p2, q2,
                                   T(ii,ii-1), T(ii-1,ii),
                                   p3, q3,
                                   T(ii,ii), tmp[0], tmp[1]
                                   ); 
        if (finish)
            break; 

        P.col(ii+1) = p3; 
        Q.col(ii+1) = q3; 
        T(ii+1,ii)= tmp[0];
        T(ii,ii+1)= tmp[1]; 
    }
    PRINT_MAT(A); 
    PRINT_MAT(P); 
    PRINT_MAT(Q); 
    PRINT_MAT(T); 
    std::cout << "P^T A Q should approximate T\n"; 
    PRINT_MAT(P.transpose()*A*Q); 
}
//##############################################################################
int main() {
    Basics();
    Tridiagonalization();
}
//##############################################################################
