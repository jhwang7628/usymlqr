#ifndef USYM_LINEAR_SOLVER_HPP
#define USYM_LINEAR_SOLVER_HPP
#include "macros.h"
#include "sparse_matrix.hpp"
#include "usym_tridiag.hpp"
//##############################################################################
// Class USYM_Linear_Solver
//##############################################################################
template<typename T, int N>
class USYM_Linear_Solver
{
    using USYM_Vector       = typename Sparse_Matrix<T,N>::USYM_Vector;
    using USYM_VectorX      = typename Sparse_Matrix<T,N>::USYM_VectorX;
    using USYM_Matrix       = typename Sparse_Matrix<T,N>::USYM_Matrix; 
    using USYM_MatrixX      = typename Sparse_Matrix<T,N>::USYM_MatrixX; 
    using Sparse_Matrix_Ptr = typename USYM_Tridiag<T,N>::Sparse_Matrix_Ptr; 

    Sparse_Matrix_Ptr _A; 
    USYM_Vector       _b; 
    std::unique_ptr<USYM_Tridiag<T,N>>  _tridiagonalization;

    // actual storage for the solver: 
    // 4 generalized Lanczos vectors, 
    // 2 one for x and one for w 
    int _lanczos_rewrite_pointer = 0;
    int _step = 0; 
    USYM_Vector _p[2]; 
    USYM_Vector _q[2];
    USYM_Vector _x; 
    USYM_Vector _w;
    USYM_Vector _betas;  
    USYM_Vector _gammas; 
    USYM_Vector _alphas; 

    // aux vecors, not strictly necessary but for convenience
    // notation following Saunders Eq 5.1 for matrix L
    USYM_Vector _deltas; 
    USYM_Vector _lambdas; 
    USYM_Vector _epsilons; 

    bool _initialized = false; 
public: 
    USYM_Linear_Solver() = default; 
    USYM_Linear_Solver(std::shared_ptr<Sparse_Matrix<T,N>> A,
                       const USYM_Vector &b)
        : _A(A),
          _b(b),
          _tridiagonalization(new USYM_Tridiag<T,N>(_A))
    {}
    void Initialize(const USYM_Vector &x0);
    void Step();
    void ComputePlaneRotation(const T &a_11, const T &a_12,
                              T &c_11, T &c_12,
                              T &c_21, T &c_22); 
    void ApplyPlaneRotation(const T &c_11, const T &c_12, 
                            const T &c_21, const T &c_22, 
                                  T &a_11,       T &a_12, 
                                  T &a_21,       T &a_22, 
                                  T &a_31,       T &a_32); // a_31 === 0 @ input
}; 

//##############################################################################
// Function Initialize
//##############################################################################
template<typename T, int N> 
void USYM_Linear_Solver<T,N>::
Initialize(const USYM_Vector &x0)
{
    assert(_A); 

    // p0, q0
    USYM_Vector &pOld = _p[0]; 
    USYM_Vector &qOld = _q[0]; 
    pOld = USYM_Vector::Zero(); 
    qOld = USYM_Vector::Zero(); 

    // p1, q1
    USYM_Vector &pNow = _p[1]; 
    USYM_Vector &qNow = _q[1]; 

    pNow = (_b - (*_A)*x0); 
    _betas(0) = pNow.norm(); 
    pNow /= _betas(0);

    qNow = USYM_Vector::Random(); 
    _gammas(0) = qNow.norm(); 
    qNow /= _gammas(0); 

    _lanczos_rewrite_pointer = 0; 
    _step = 0;
    _initialized = true; 
}

//##############################################################################
// Function Step
//##############################################################################
template<typename T, int N> 
void USYM_Linear_Solver<T,N>::
Step()
{
    assert(_initialized && _tridiagonalization); 
    if (_step >= N) // silently return
        return; 
    USYM_Vector &pOld = _p[ _lanczos_rewrite_pointer     ]; 
    USYM_Vector &qOld = _q[ _lanczos_rewrite_pointer     ]; 
    USYM_Vector &pNow = _p[(_lanczos_rewrite_pointer+1)%2]; 
    USYM_Vector &qNow = _q[(_lanczos_rewrite_pointer+1)%2]; 

    T tmp[2];  // beta, gamma
    bool done = true; 
    if (_step > 0)
        done = _tridiagonalization->StepInPlace(pOld, qOld,
                                                pNow, qNow, 
                                                _betas[_step], _gammas[_step], 
                                                _alphas[_step], tmp[0], tmp[1]); 
    if (!done && _step<N-1)
    {
        _betas[_step+1]  = tmp[0]; 
        _gammas[_step+1] = tmp[1];
    }
    else
    {
        done = true; 
    }

    // for step i, do plane rotation on step i-1
    _deltas[_step] = _alphas[_step]; 
    if (!done && _step>0)
    {
        _lambdas[_step] = _betas[_step+1]; 
        _epsilons[_step] = 0.; 
        T c[4]; // c_11, c_12, c_21, c_22
        ComputePlaneRotation(_deltas[_step-1], _gammas[_step],
                             c[0], c[1], 
                             c[2], c[3]); 
        T upper_tmp = _gammas[_step]; 
        ApplyPlaneRotation(c[0], c[1],
                           c[2], c[3], 
                           _deltas[_step-1] , upper_tmp, 
                           _lambdas[_step-1], _deltas[_step], 
                           _epsilons[_step-1], _lambdas[_step]); 
        PRINT(std::cout, upper_tmp); // should be close to zero
    }

    _lanczos_rewrite_pointer = (_lanczos_rewrite_pointer+1)%2; 
    ++_step; 
}

//##############################################################################
// Function ComputePlaneRotation
//##############################################################################
template<typename T, int N> 
void USYM_Linear_Solver<T,N>::
ComputePlaneRotation(const T &a_11, const T &a_12,
                     T &c_11, T &c_12,
                     T &c_21, T &c_22)
{
    // this means tan(theta) = c_12/c_22 = -a_12/a_11
    const T theta = atan(-a_12/a_11); 
    c_11 = cos(theta); c_12 = sin(theta); 
    c_21 =-sin(theta); c_22 = cos(theta); 
}
//##############################################################################
// Function ApplyPlaneRotation
//
// [[a_11, a_12],     [[c_11, c_12],
//  [a_21, a_22],  *   [c_21, c_22]]
//  [a_31, a_32]]
//##############################################################################
template<typename T, int N> 
void USYM_Linear_Solver<T,N>::
ApplyPlaneRotation(const T &c_11, const T &c_12, 
                   const T &c_21, const T &c_22, 
                         T &a_11,       T &a_12, 
                         T &a_21,       T &a_22, 
                         T &a_31,       T &a_32) // a_31 === 0 @ input
{
    a_31 = (T)0; //invariant
    T a_11_n = a_11*c_11 + a_12*c_21; 
    T a_12_n = a_11*c_12 + a_12*c_22; 
    T a_21_n = a_21*c_11 + a_22*c_21; 
    T a_22_n = a_21*c_12 + a_22*c_22; 
    T a_31_n = a_31*c_11 + a_32*c_21; 
    T a_32_n = a_31*c_12 + a_32*c_22; 
    std::swap(a_11, a_11_n); 
    std::swap(a_12, a_12_n); 
    std::swap(a_21, a_21_n); 
    std::swap(a_22, a_22_n); 
    std::swap(a_31, a_31_n); 
    std::swap(a_32, a_32_n); 
}
#endif
