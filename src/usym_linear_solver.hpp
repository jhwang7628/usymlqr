#ifndef USYM_LINEAR_SOLVER_HPP
#define USYM_LINEAR_SOLVER_HPP
#include "macros.h"
#include "sparse_matrix.hpp"
#include "usym_tridiag.hpp"
#include "tridiagonal_matrix.hpp"
#include "lower_triangular_matrix.hpp"
//##############################################################################
// Class USYM_Linear_Solver
//##############################################################################
template<typename T, class T_Vector, class T_Matrix>
class USYM_Linear_Solver
{
    std::shared_ptr<T_Matrix> _A; 
    T_Vector       _b; 
    std::unique_ptr<USYM_Tridiag<T,T_Vector,T_Matrix>>  _tridiagonalization;
    int _N;

    // actual storage for the solver: 
    // 4 generalized Lanczos vectors, 
    // 2 one for x and one for w 
    int _lanczos_rewrite_pointer = 0;
    int _step = 0; 
    T_Vector _p[2]; 
    T_Vector _q[2];
    T_Vector _x; 
    T_Vector _w;
    Tridiagonal_Matrix<T> _matT; 
    Lower_Triangular_Matrix<T> _matL; 
    T_Vector _z; 

    bool _initialized = false; 
public: 
    USYM_Linear_Solver() = default; 
    USYM_Linear_Solver(std::shared_ptr<T_Matrix> A,
                       const T_Vector &b)
        : _A(A),
          _b(b),
          _tridiagonalization(new USYM_Tridiag<T,T_Vector,T_Matrix>(_A)),
          _N(b.size()),
          _p{T_Vector(_N),T_Vector(_N)},
          _q{T_Vector(_N),T_Vector(_N)},
          _x(T_Vector(_N)),
          _w(T_Vector(_N)),
          _matT(_N),
          _matL(_N),
          _z(T_Vector(_N))
    {}
    T_Vector Initialize(const T_Vector &x0);
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
template<typename T, class T_Vector, class T_Matrix>
T_Vector USYM_Linear_Solver<T,T_Vector,T_Matrix>::
Initialize(const T_Vector &x0)
{
    assert(_A); 

    // p0, q0
    T_Vector &pOld = _p[0]; 
    T_Vector &qOld = _q[0]; 
    pOld.setZero(); 
    qOld.setZero(); 

    // p1, q1
    T_Vector &pNow = _p[1]; 
    T_Vector &qNow = _q[1]; 

    T beta_1, gamma_1; 
    T_Vector r0 = (_b - (*_A)*x0);
    beta_1 = r0.norm(); 
    pNow = r0 / beta_1;

    qNow = T_Vector::Random(_N); 
    gamma_1 = qNow.norm(); 
    qNow /= gamma_1; 

    int n = _matT.AddBetaAndGamma(beta_1,gamma_1); 
    assert(n==1); 

    // initialize x, w
    _x = x0; 
    _w.setZero(); 
    _z.setZero(); 

    _lanczos_rewrite_pointer = 0; 
    _step = 0;
    _initialized = true; 

    return r0; 
}

//##############################################################################
// Function Step
//##############################################################################
template<typename T, class T_Vector, class T_Matrix>
void USYM_Linear_Solver<T,T_Vector,T_Matrix>::
Step()
{
    const int &i = _step; 
    std::cout << "========== STEP " << i << " ========== " << std::endl;
    assert(_initialized && _tridiagonalization); 
    T_Vector &pOld = _p[ _lanczos_rewrite_pointer     ]; 
    T_Vector &qOld = _q[ _lanczos_rewrite_pointer     ]; 
    T_Vector &pNow = _p[(_lanczos_rewrite_pointer+1)%2]; 
    T_Vector &qNow = _q[(_lanczos_rewrite_pointer+1)%2]; 

    std::cout << " step tridiagonalization\n";
    const T beta_i  = _matT(i  , i-1); 
    const T gamma_i = _matT(i-1, i  ); 
    T alpha_i; // to be determined
    T beta_ip1, gamma_ip1;  // beta, gamma
    bool done = _tridiagonalization->StepInPlace(pOld, qOld,
                                                 pNow, qNow, 
                                                 beta_i, gamma_i, 
                                                 alpha_i, beta_ip1, gamma_ip1); 
    {
        const int n = _matT.AddAlpha(alpha_i); 
        assert(n == i+1); 
    }
    if (!done)
    {
        const int n = _matT.AddBetaAndGamma(beta_ip1, gamma_ip1); 
        assert(n == i+2); 
    }

    // initialize matrix L in first two steps
    if (i == 1) 
    {
        _matL.AddColumn((T)0., _matT(0,0), _matT(1,0), _matT(2,0));
        _matL.AddColumn(_matT(0,1), _matT(1,1), _matT(2,1), (T)0.); 
    } 
    else if (i > 1 && !done)
    {
        _matL.AddColumn(_matT(i-1,i), _matT(i,i), _matT(i+1,i), (T)0.); 
    }

    // for step i, do plane rotation on (i-1,i-1), (i-1,i), (i,i-1), (i,i)
    // so that element (i-1,i)=0
    if (!done && i>0)
    {
        T &a11 = _matL(i-1,i-1); 
        T &a12 = _matL(i-1,i  ); 
        T &a21 = _matL(i  ,i-1); 
        T &a22 = _matL(i  ,i  ); 
        T  a31 = 0.0;  // software design..
        T &a32 = _matL(i+1,i  ); 
        assert(fabs(a31-0)<SMALL_NUM); // a31 should be zero before rotation

        T c11, c12, c21, c22; 
        std::cout << "T = \n";
        _matT.Print(i+2, i+1); 
        std::cout << "L = \n";
        _matL.Print(i+2, i+1);
        ComputePlaneRotation(a11, a12,
                             c11, c12, 
                             c21, c22); 
        std::cout << " plane rotation\n";
        ApplyPlaneRotation(c11, c12,
                           c21, c22, 
                           a11, a12,
                           a21, a22,
                           a31, a32); 
        assert(fabs(a12-0)<SMALL_NUM); // a12 should be zero after rotation
        _matL(i+1,i-1) = a31;  // write back
        std::cout << "L = \n";
        _matL.Print(i+2, i+1);
    }

    // update x and w
    //const T_Vector &q_i = qNow; // just to make life easier with alias..
    //if (_step==0)
    //    _z(_step) = _betas(_step) / _deltas(_step); 
    //else if (_step==1)
    //    _z(_step) = - _lambdas(_step-1)*_z(_step-1) / _deltas(_step); 
    //else 
    //    _z(_step) = -(_lambdas(_step-1)*_z(_step-1) 
    //                + _epsilons(_step-2)*_z(_step-2)) / _deltas(_step); 

    // TODO here

    _lanczos_rewrite_pointer = (_lanczos_rewrite_pointer+1)%2; 
    ++_step; 
}

//##############################################################################
// Function ComputePlaneRotation
//##############################################################################
template<typename T, class T_Vector, class T_Matrix>
void USYM_Linear_Solver<T,T_Vector,T_Matrix>::
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
template<typename T, class T_Vector, class T_Matrix>
void USYM_Linear_Solver<T,T_Vector,T_Matrix>::
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
