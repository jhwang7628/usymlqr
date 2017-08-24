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
    T_Vector                  _b; 
    std::unique_ptr<USYM_Tridiag<T,T_Vector,T_Matrix>>  _tridiagonalization;
    int _N;

    // solver config
    int _maxItn;
    T _a_tol; 
    T _b_tol; 

    // actual storage for the solver: 
    // 4 generalized Lanczos vectors, 
    // 2 one for x and one for w 
    int _lanczos_rewrite_pointer = 0;
    int _step = 0; 
    T_Vector _p[3]; 
    T_Vector _q[3];
    T_Vector _x; 
    T_Vector _w;
    Tridiagonal_Matrix<T> _matT; 
    Lower_Triangular_Matrix<T> _matL; 

    // cache
    T _bnorm; 
    T _rhs_1;  // one step ago
    T _rhs_2;  // two step ago
    T _Anorm2; // norm(A) ~= norm(T)
    T _Anorm; 
    T _xnorm2; // norm(x)  = norm(h_j) when x0 = 0
    T _xnorm; 
    T _zbar; 
    T _rnorm; 
    T _eta; 
    T _cgnorm;
    std::vector<T> _z;  // need fast dynamic push_back and memory manage

    bool _initialized = false; 

public: 
    enum Mode{USYMLQ=0} mode = USYMLQ; 
    
    USYM_Linear_Solver() = default; 
    USYM_Linear_Solver(std::shared_ptr<T_Matrix> A,
                       const T_Vector &b)
        : _A(A),
          _b(b),
          _tridiagonalization(new USYM_Tridiag<T,T_Vector,T_Matrix>(_A)),
          _N(b.size()),
          _maxItn(_N),
          _a_tol(DEFAULT_TOL),
          _b_tol(DEFAULT_TOL),
          _p{T_Vector(_N),T_Vector(_N)},
          _q{T_Vector(_N),T_Vector(_N)},
          _x(T_Vector(_N)),
          _w(T_Vector(_N)),
          _matT(_N),
          _matL(_N)
    {
        _z.reserve(_N); 
    }
    inline void SetMaxIteration(const int maxN) 
        {_maxItn = maxN;}
    inline void SetTol(const T atol, const T btol)
        {_a_tol = atol; _b_tol = btol;}
    T_Vector Initialize(const T_Vector &x0);
    int Solve(T_Vector &x, T &rnorm); 

private: 
    int Step();
    void ComputePlaneRotation(const T &a_11, const T &a_12,
                              T &c_11, T &c_12,
                              T &c_21, T &c_22); 
    void ApplyPlaneRotation(const T &c_11, const T &c_12, 
                            const T &c_21, const T &c_22, 
                                  T &a_11,       T &a_12, 
                                  T &a_21,       T &a_22, 
                                  T &a_31,       T &a_32); // a_31 === 0 @ input
    //// debug /////
    void CheckError_z(); 
    void CheckResidual(); 
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
    _w = qNow; 
    _rhs_1 = beta_1; 
    _rhs_2 = 0.0;
    _Anorm2 = 0.0; 
    _xnorm2 = 0.0; 
    _bnorm = _b.norm();

    _lanczos_rewrite_pointer = 2; 
    _step = 0;
    _initialized = true; 

    return r0; 
}

//##############################################################################
// Function Solve
//##############################################################################
template<typename T, class T_Vector, class T_Matrix>
int USYM_Linear_Solver<T,T_Vector,T_Matrix>::
Solve(T_Vector &x, T&rnorm)
{
    assert(_initialized); 
    int flag = 0;
    while (flag == 0)
        flag = Step(); 
    x = _x; 
    rnorm = _rnorm; 
    return flag; 
}

//##############################################################################
// Function Step
//##############################################################################
template<typename T, class T_Vector, class T_Matrix>
int USYM_Linear_Solver<T,T_Vector,T_Matrix>::
Step()
{
    const int &i = _step; 
    int flag = 0;
    std::cout << "STEP " << i;
    assert(_initialized && _tridiagonalization); 
    T_Vector &p_im1 = _p[(_lanczos_rewrite_pointer+1)%3]; 
    T_Vector &q_im1 = _q[(_lanczos_rewrite_pointer+1)%3]; 
    T_Vector &p_i   = _p[(_lanczos_rewrite_pointer+2)%3]; 
    T_Vector &q_i   = _q[(_lanczos_rewrite_pointer+2)%3]; 
    T_Vector &p_ip1 = _p[ _lanczos_rewrite_pointer]; 
    T_Vector &q_ip1 = _q[ _lanczos_rewrite_pointer]; 

    const T beta_i  = _matT(i  , i-1); 
    const T gamma_i = _matT(i-1, i  ); 
    T alpha_i; // to be determined
    T beta_ip1, gamma_ip1;  // beta, gamma
    _tridiagonalization->Step(p_im1, q_im1,
                                          p_i  , q_i  , 
                                          beta_i, gamma_i, 
                                          p_ip1, q_ip1,
                                          alpha_i, beta_ip1, gamma_ip1); 
    int ntest;
    ntest = _matT.AddAlpha(alpha_i); assert(ntest == i+1); 
    ntest = _matT.AddBetaAndGamma(beta_ip1, gamma_ip1); assert(ntest == i+2); 

    // initialize matrix L in first two steps
    if (i == 1) 
    {
        _matL.AddColumn((T)0., _matT(0,0), _matT(1,0), _matT(2,0));
        _matL.AddColumn(_matT(0,1), _matT(1,1), _matT(2,1), (T)0.); 
    } 
    else if (i > 1)
    {
        _matL.AddColumn(_matT(i-1,i), _matT(i,i), _matT(i+1,i), (T)0.); 
    }

    // for step i, do plane rotation on (i-1,i-1), (i-1,i), (i,i-1), (i,i)
    // so that element (i-1,i)=0
    if (i > 0)
    {
        T &a11 = _matL(i-1,i-1); 
        T &a12 = _matL(i-1,i  ); 
        T &a21 = _matL(i  ,i-1); 
        T &a22 = _matL(i  ,i  ); 
        T  a31 = 0.0;  // software design..
        T &a32 = _matL(i+1,i  ); 
        assert(fabs(a31-0)<SMALL_NUM); // a31 should be zero before rotation

        T c11, c12, c21, c22; 
        ComputePlaneRotation(a11, a12,
                             c11, c12, 
                             c21, c22); 
        ApplyPlaneRotation(c11, c12,
                           c21, c22, 
                           a11, a12,
                           a21, a22,
                           a31, a32); 
        assert(fabs(a12-0)<SMALL_NUM); // a12 should be zero after rotation
        _matL(i+1,i-1) = a31;  // write back

        const T z = _rhs_1 / _matL(i-1,i-1); 
        const T s = z*c11; 
        const T t = z*c21; 

        _x += _w*s   + q_i*t  ; 
        _w  = _w*c12 + q_i*c22; 
        _rhs_1 = _rhs_2 - _matL(i  , i-1)*z; 
        _rhs_2 =        - _matL(i+1, i-1)*z; 
        _z.push_back(z); 

        _zbar = _rhs_1 / _matL(i,i); // estimate of z_i
        _eta  = c21*z + c22*_zbar;    // last element of h_j
        _cgnorm = _matT(i+1,i) * fabs(_eta); 
        _Anorm2 += pow(_matT(i  ,i  ),2) 
                 + pow(_matT(i  ,i+1),2)
                 + pow(_matT(i+1,i  ),2); 
        _Anorm = sqrt(_Anorm2); 
        _xnorm = sqrt(_xnorm2 + pow(_zbar,2)); 
        _xnorm2 += pow(z,2); 

        const T t_rel = _cgnorm / _bnorm; 
        const T t_tol = _b_tol + _a_tol*_Anorm*_xnorm/_bnorm; 
        const T t1    = 1.0 + t_rel / (1.0 + _Anorm*_xnorm/_bnorm);
        if      (t_rel < t_tol) flag = 1; 
        else if (t1 <= 1.0)     flag = 4; 

        //CheckError_z();
        //CheckResidual(); 
    }
    std::cout << "   norm(r_cg) = " << _cgnorm << std::endl;
    _lanczos_rewrite_pointer = (_lanczos_rewrite_pointer+1)%3; 
    if (_step++ >= _maxItn) flag = 7; 
    if (flag != 0) 
    {
        _x = _x + _zbar*_w; 
        _rnorm = _cgnorm; 
    }
    return flag; 
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

//##############################################################################
// Function CheckError_z
//##############################################################################
template<typename T, class T_Vector, class T_Matrix>
void USYM_Linear_Solver<T,T_Vector,T_Matrix>::
CheckError_z()
{
    const int N = _z.size(); 
    std::vector<T> beta(N, 0.0); 
    for (int r=0; r<N; ++r)
    for (int c=0; c<N; ++c)
        beta.at(r) += _z.at(c)*_matL(r,c); 
    beta.at(0) -= _matT.Get_Beta_1(); 
    std::cout << "error for solving z = ";
    std::copy(beta.begin(),beta.end(),std::ostream_iterator<T>(std::cout," ")); 
    std::cout << std::endl;
}

//##############################################################################
// Function CheckResidual
//##############################################################################
template<typename T, class T_Vector, class T_Matrix>
void USYM_Linear_Solver<T,T_Vector,T_Matrix>::
CheckResidual()
{
    T_Vector r = (*_A)*_x - _b; 
    const int N = _b.size(); 
    const T_Vector xstar = _A->fullPivLu().solve(_b); 
    std::cout << "X^*  = ";
    for (int ii=0; ii<N; ++ii)
        printf("% 8.4f ", xstar[ii]); 
    std::cout << std::endl; 
    std::cout << "X_cg = ";
    for (int ii=0; ii<N; ++ii)
        printf("% 8.4f ", _x[ii]); 
    std::cout << std::endl; 
    std::cout << "Res  = ";
    for (int ii=0; ii<N; ++ii)
        printf("% 8.4f ", r[ii]); 
    std::cout << std::endl; 
    std::cout << "||Res|| = " << r.norm() << std::endl;
}
#endif
