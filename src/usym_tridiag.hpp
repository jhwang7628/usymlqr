#ifndef USYM_TRIDIAG_HPP
#define USYM_TRIDIAG_HPP
#include <memory>
#include "macros.h"
#include "sparse_matrix.hpp"
//##############################################################################
// Class USym_Tridiag
//   A class for unsymmetric tridiagonalization, described in Saunders et al.
//   "Two Conjugate-gradient-type Methods for Unsymmetric Linear Equations"
//   Algorithm 1. 
//##############################################################################
template<typename T, int N>
class USYM_Tridiag
{
public: 
    using USYM_Vector = typename Sparse_Matrix<T,N>::USYM_Vector;
    using USYM_VectorX= typename Sparse_Matrix<T,N>::USYM_VectorX;
    using USYM_Matrix = typename Sparse_Matrix<T,N>::USYM_Matrix; 
    using USYM_MatrixX= typename Sparse_Matrix<T,N>::USYM_MatrixX; 
    using Sparse_Matrix_Ptr = std::shared_ptr<Sparse_Matrix<T,N>>; 

private: 
    Sparse_Matrix_Ptr _A; 
    // buffer
    USYM_Vector _u; 
    USYM_Vector _v; 

public:
    USYM_Tridiag() = default; 
    USYM_Tridiag(const Sparse_Matrix_Ptr &A)
        : _A(A)
    {}
    void Set_A(const Sparse_Matrix_Ptr &A){_A = A;}
    // i=1 step needs special initialization
    bool InitialStep(const USYM_Vector &b, const USYM_Vector &c,
                     USYM_Vector &p_1, USYM_Vector &q_1, // output
                     USYM_Vector &p_2, USYM_Vector &q_2, // output
                     T &alpha_1, T &beta_2, T &gamma_2
                    );
    bool Step(const USYM_Vector &p_im1, const USYM_Vector &q_im1,
              const USYM_Vector &p_i  , const USYM_Vector &q_i  ,
              const T beta_i, const T gamma_i,
              USYM_Vector &p_ip1, USYM_Vector &q_ip1, // output
              T &alpha_i, T &beta_ip1, T &gamma_ip1); 
    bool StepInPlace(USYM_Vector &p_im1, USYM_Vector &q_im1, // rewrite im1 
                     const USYM_Vector &p_i  , const USYM_Vector &q_i  ,
                     const T beta_i, const T gamma_i,
                     T &alpha_i, T &beta_ip1, T &gamma_ip1); 
};

//##############################################################################
// Function InitialStep
//##############################################################################
template<typename T, int N> 
bool USYM_Tridiag<T,N>::
InitialStep(const USYM_Vector &b, const USYM_Vector &c,
            USYM_Vector &p_1, USYM_Vector &q_1, // output
            USYM_Vector &p_2, USYM_Vector &q_2,
            T &alpha_1, T &beta_2, T &gamma_2
           )
{
    USYM_Vector p_0 = USYM_Vector::Zero(); 
    USYM_Vector q_0 = USYM_Vector::Zero(); 
    const T beta_1  = b.norm(); 
    const T gamma_1 = c.norm(); 
    p_1 = b/beta_1; 
    q_1 = c/gamma_1;
    // i=1 step
    return Step(p_0, q_0,
                p_1, q_1, 
                beta_1, gamma_1,
                p_2, q_2, 
                alpha_1, beta_2, gamma_2); 
}

//##############################################################################
// Function Step
//##############################################################################
template<typename T, int N> 
bool USYM_Tridiag<T,N>::
Step(const USYM_Vector &p_im1, const USYM_Vector &q_im1,
     const USYM_Vector &p_i  , const USYM_Vector &q_i ,
     const T beta_i, const T gamma_i,
     USYM_Vector &p_ip1, USYM_Vector &q_ip1, // output
     T &alpha_i, T &beta_ip1, T &gamma_ip1)
{
    assert(_A); 
    // u = A q_i - gamma_i p_{i-1}
    _A->Premultiply_By_Matrix(q_i, _u); 
    _u -= (gamma_i * p_im1); 
    // v = A^* p_i - beta_i q_{i-1}
    _A->Premultiply_By_Matrix_Conjugate(p_i, _v); 
    _v -= (beta_i * q_im1); 
    // alpha = p_i^T u
    alpha_i = p_i.transpose()*_u; 
    // continue..
    _u -= (alpha_i * p_i); 
    _v -= (alpha_i * q_i); 
    beta_ip1  = _u.norm(); 
    gamma_ip1 = _v.norm(); 
    if (fabs(beta_ip1) < 1E-8 || fabs(gamma_ip1) < 1E-8)
        return true; 
    p_ip1 = _u / beta_ip1;
    q_ip1 = _v / gamma_ip1; 
    return false; 
}

//##############################################################################
// Function StepInPlace
//##############################################################################
template<typename T, int N> 
bool USYM_Tridiag<T,N>::
StepInPlace(USYM_Vector &p_im1, USYM_Vector &q_im1,
            const USYM_Vector &p_i  , const USYM_Vector &q_i ,
            const T beta_i, const T gamma_i,
            T &alpha_i, T &beta_ip1, T &gamma_ip1)
{
    assert(_A); 
    // u = A q_i - gamma_i p_{i-1}
    _A->Premultiply_By_Matrix(q_i, _u); 
    _u -= (gamma_i * p_im1); 
    // v = A^* p_i - beta_i q_{i-1}
    _A->Premultiply_By_Matrix_Conjugate(p_i, _v); 
    _v -= (beta_i * q_im1); 
    // alpha = p_i^T u
    alpha_i = p_i.transpose()*_u; 
    // continue..
    _u -= (alpha_i * p_i); 
    _v -= (alpha_i * q_i); 
    beta_ip1  = _u.norm(); 
    gamma_ip1 = _v.norm(); 
    if (fabs(beta_ip1) < 1E-8 || fabs(gamma_ip1) < 1E-8)
        return true; 
    p_im1 = _u / beta_ip1;
    q_im1 = _v / gamma_ip1; 
    return false; 
}

#endif
