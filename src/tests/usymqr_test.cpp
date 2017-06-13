#include <memory>
#include <iostream>
#include "sparse_matrix.hpp"
#include "usym_tridiag.hpp"
#include "Eigen/Dense"
//##############################################################################
void Basics() {
    Eigen::Matrix3d m;
    m << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    std::cout << m << std::endl;

    Sparse_Matrix<double,10> sm; 
}
//##############################################################################
void Tridiagonalization() {
    const int N = 10; 
    using USYM_Vector = typename Sparse_Matrix<double,N>::USYM_Vector;
    using USYM_VectorX= typename Sparse_Matrix<double,N>::USYM_VectorX;
    using USYM_Matrix = typename Sparse_Matrix<double,N>::USYM_Matrix; 
    using USYM_MatrixX= typename Sparse_Matrix<double,N>::USYM_MatrixX; 

    USYM_Matrix m = USYM_Matrix::Random();
    auto sm = std::make_unique<Sparse_Matrix<double,10>>(m); 
    USym_Tridiag<double,N> tridiag(sm); 
    auto b = USYM_Vector::Random(N); 
    auto c = USYM_Vector::Random(N); 
    tridiag.Initialize(b, c);
}
//##############################################################################
int main() {
    Basics();
    Tridiagonalization();
}
//##############################################################################
