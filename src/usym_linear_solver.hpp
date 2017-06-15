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
    T _initial_beta; 
    T _initial_gamma; 
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
    _initial_beta = pNow.norm(); 
    pNow /= _initial_beta;

    qNow = USYM_Vector::Random(); 
    _initial_gamma = qNow.norm(); 
    qNow /= _initial_gamma; 

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

    T tmp[2]; 
    bool done; 
    if (_step > 0)
        done = _tridiagonalization->StepInPlace(pOld, qOld,
                                                pNow, qNow, 
                                                _betas[_step], _gammas[_step], 
                                                _alphas[_step], tmp[0], tmp[1]); 
    else
        done = _tridiagonalization->StepInPlace(pOld, qOld,
                                                pNow, qNow, 
                                                _initial_beta, _initial_gamma, 
                                                _alphas[_step], tmp[0], tmp[1]); 
    if (!done && _step<N-1)
    {
        _betas[_step+1]  = tmp[0]; 
        _gammas[_step+1] = tmp[1];
    }
    else
    {
        std::cout << "done!\n";
        return; 
    }

    // FIXME debug START
    std::cout << "+++ Step " << _step << ": \n"; 
    std::cout << "a,b,g = " << _alphas[_step] << ", " 
                            << tmp[0]         << ", "
                            << tmp[1]         << std::endl;
    // FIXME debug END
    _lanczos_rewrite_pointer = (_lanczos_rewrite_pointer+1)%2; 
    ++_step; 
}
#endif
