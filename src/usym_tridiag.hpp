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
template<typename T, class T_Vector, class T_Matrix>
class USYM_Tridiag
{
private: 
    std::shared_ptr<T_Matrix> _A; 
    // buffer
    T_Vector _u; 
    T_Vector _v; 
    // initial vectors
    T_Vector _b; 
    T_Vector _c; 

public:
    USYM_Tridiag() = default; 
    USYM_Tridiag(const std::shared_ptr<T_Matrix> &A)
        : _A(A), _u(A->rows()), _v(A->rows())
    {
        _b = T_Vector::Random(A->rows()); 
        _c = T_Vector::Random(A->rows()); 
    }
    void Set_A(const std::shared_ptr<T_Matrix> &A){_A = A;}
    void Set_b(const T_Vector &b){_b = b;}
    void Set_c(const T_Vector &c){_c = c;}
    // i=1 step needs special initialization
    bool InitialStep(T_Vector &p_1, T_Vector &q_1, // output
                     T_Vector &p_2, T_Vector &q_2, // output
                     T &alpha_1, T &beta_2, T &gamma_2
                    );
    bool Step(const T_Vector &p_im1, const T_Vector &q_im1,
              const T_Vector &p_i  , const T_Vector &q_i  ,
              const T beta_i, const T gamma_i,
              T_Vector &p_ip1, T_Vector &q_ip1, // output
              T &alpha_i, T &beta_ip1, T &gamma_ip1); 
    bool StepInPlace(T_Vector &p_im1, T_Vector &q_im1, // rewrite im1 
                     const T_Vector &p_i  , const T_Vector &q_i  ,
                     const T beta_i, const T gamma_i,
                     T &alpha_i, T &beta_ip1, T &gamma_ip1); 
};

//##############################################################################
// Function InitialStep
//##############################################################################
template<typename T, class T_Vector, class T_Matrix>
bool USYM_Tridiag<T,T_Vector,T_Matrix>::
InitialStep(T_Vector &p_1, T_Vector &q_1, // output
            T_Vector &p_2, T_Vector &q_2,
            T &alpha_1, T &beta_2, T &gamma_2
           )
{
    T_Vector p_0 = T_Vector::Zero(_A->rows()); 
    T_Vector q_0 = T_Vector::Zero(_A->rows()); 
    const T beta_1  = _b.norm(); 
    const T gamma_1 = _c.norm(); 
    p_1 = _b/beta_1; 
    q_1 = _c/gamma_1;
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
template<typename T, class T_Vector, class T_Matrix>
bool USYM_Tridiag<T,T_Vector,T_Matrix>::
Step(const T_Vector &p_im1, const T_Vector &q_im1,
     const T_Vector &p_i  , const T_Vector &q_i ,
     const T beta_i, const T gamma_i,
     T_Vector &p_ip1, T_Vector &q_ip1, // output
     T &alpha_i, T &beta_ip1, T &gamma_ip1)
{
    assert(_A); 
    // u = A q_i - gamma_i p_{i-1}
    _u = (*_A)*q_i - gamma_i*p_im1; 
    // v = A^* p_i - beta_i q_{i-1}
    _v = _A->transpose()*p_i - beta_i*q_im1; 
    // alpha = p_i^T u
    alpha_i = p_i.dot(_u); 
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
template<typename T, class T_Vector, class T_Matrix>
bool USYM_Tridiag<T,T_Vector,T_Matrix>::
StepInPlace(T_Vector &p_im1, T_Vector &q_im1,
            const T_Vector &p_i  , const T_Vector &q_i ,
            const T beta_i, const T gamma_i,
            T &alpha_i, T &beta_ip1, T &gamma_ip1)
{
    assert(_A); 
    // u = A q_i - gamma_i p_{i-1}
    _u = (*_A)*q_i - gamma_i*p_im1; 
    // v = A^* p_i - beta_i q_{i-1}
    _v = _A->transpose()*p_i - beta_i*q_im1; 
    // alpha = p_i^T u
    alpha_i = p_i.dot(_u); 
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
