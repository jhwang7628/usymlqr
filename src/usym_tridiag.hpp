#ifndef USYM_TRIDIAG_HPP
#define USYM_TRIDIAG_HPP
#include <memory>
#include "sparse_matrix.hpp"
//##############################################################################
// Class USym_Tridiag
//   A class for unsymmetric tridiagonalization, described in Saunders et al.
//   "Two Conjugate-gradient-type Methods for Unsymmetric Linear Equations"
//   Algorithm 1. 
//##############################################################################
template<typename T, int N>
class USym_Tridiag
{
    using USYM_Vector = typename Sparse_Matrix<T,N>::USYM_Vector;
    using USYM_VectorX= typename Sparse_Matrix<T,N>::USYM_VectorX;
    using USYM_Matrix = typename Sparse_Matrix<T,N>::USYM_Matrix; 
    using USYM_MatrixX= typename Sparse_Matrix<T,N>::USYM_MatrixX; 
    using Sparse_Matrix_Ptr = std::unique_ptr<Sparse_Matrix<T,N>>; 
  
    Sparse_Matrix_Ptr _A; 
    USYM_MatrixX _Q; 
    USYM_MatrixX _P; 
    USYM_VectorX _betas;
    USYM_VectorX _gammas; 
    int _i = 0; 
public:
    USym_Tridiag() = default; 
    USym_Tridiag(Sparse_Matrix_Ptr &A)
        : _A(std::move(A))
    {}
    void Set_A(const Sparse_Matrix_Ptr &A){_A = std::move(A);}
    void Initialize(const USYM_Vector &b, const USYM_Vector &c);
};

//##############################################################################
// Function Initialize
//##############################################################################
template<typename T, int N> 
void USym_Tridiag<T,N>::
Initialize(const USYM_Vector &b, const USYM_Vector &c)
{
    assert(_A); 
    const USYM_Vector p_m1 = USYM_Vector::Zero(N,1); 
    const USYM_Vector q_m1 = USYM_Vector::Zero(N,1); 
    const T beta_0  = b.norm(); 
    const T gamma_0 = c.norm(); 
    _Q.resize(Eigen::NoChange, 1);
    _P.resize(Eigen::NoChange, 1); 
    _Q.col(0) = b/beta_0; 
    _P.col(0) = c/gamma_0; 
}

#endif
