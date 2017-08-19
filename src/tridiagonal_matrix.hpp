#ifndef TRIDIAGONAL_MATRIX_HPP
#define TRIDIAGONAL_MATRIX_HPP
#include <vector>
#include "macros.h" 
//##############################################################################
// Class Tridiagonal_Matrix
//          |-   a_1 g_2  0 ................ 0     -|
//          |    b_2 a_2 g_3  0 ............ 0      |
//          |     0  b_3 a_3 g_4  0 ........ 0      |
//   T =    |     0   0  b_4 a_4 g_5  0 .... 0      |
//          |     .           .   .  .   ... 0      |
//          |     .              .   .              |
//          |     .                  .   .  g_n     |
//          |-    0   .   .   .   .    b_n  a_n    -| 
//
//##############################################################################
template<typename T>
class Tridiagonal_Matrix
{
    std::vector<T> _betas;   // b
    std::vector<T> _gammas;  // g
    std::vector<T> _alphas;  // a
    T _ZERO = 0.0;
public: 
    Tridiagonal_Matrix(const int N_hint = 1)
    {
        assert(N_hint > 0);
        _betas.reserve(N_hint); 
        _gammas.reserve(N_hint);
        _alphas.reserve(N_hint);
    }
    inline int AddAlpha(const T &alpha)
    {
        _alphas.push_back(alpha); return _alphas.size();
    }
    inline int AddBetaAndGamma(const T &beta, const T &gamma)
    {
        _betas.push_back(beta);
        _gammas.push_back(gamma);
        assert(_betas.size() == _gammas.size()); 
        return _betas.size();
    }
    inline T &operator ()(const int row, const int col)
    {
        if (row== 0 && col==-1) return _betas.at(0);
        if (row==-1 && col== 0) return _gammas.at(0);
        assert(row>=0 && col>=0);
        const int d = row - col; 
        if      (d== 0) return _alphas.at(row); 
        else if (d== 1) return _betas.at(row); 
        else if (d==-1) return _gammas.at(col);
        else           {_ZERO=(T)0.0; return _ZERO;}
    }
    inline T operator ()(const int row, const int col) const
    {
        T val = this->operator()(row,col); 
        return val; 
    }
    ///// debug methods /////
    void Print();
    void Print(const int N);
    void Print(const int Nrow, const int Ncol);
}; 

//##############################################################################
// Function Print
//##############################################################################
template<typename T>
void Tridiagonal_Matrix<T>::
Print()
{
    std::copy(_gammas.begin(), _gammas.end(),
              std::ostream_iterator<T>(std::cout," ")); 
    std::cout << std::endl; 
    std::copy(_alphas.begin(), _alphas.end(),
              std::ostream_iterator<T>(std::cout," ")); 
    std::cout << std::endl; 
    std::copy(_betas.begin(), _betas.end(),
              std::ostream_iterator<T>(std::cout," ")); 
    std::cout << std::endl; 
}

template<typename T>
void Tridiagonal_Matrix<T>::
Print(const int N)
{
    Print(N,N); 
}

template<typename T>
void Tridiagonal_Matrix<T>::
Print(const int Nrow, const int Ncol)
{
    for (int r=0; r<Nrow; ++r)
    {
        for (int c=0; c<Ncol; ++c)
            printf("% 8.4f ", (*this)(r,c)); 
        std::cout << std::endl;
    }
}
#endif
