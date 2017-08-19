#ifndef LOWER_TRIANGULAR_MATRIX_HPP
#define LOWER_TRIANGULAR_MATRIX_HPP
#include "macros.h" 
//##############################################################################
// Class Lower_Triangular_Matrix
//          |-   d_1 s_2  0 ................   0    -|
//          |    l_1 d_2 s_3  0 ............   0     |
//          |    e_1 l_2 d_3 s_4  0 ........   0     |
//   T =    |     0  e_2 l_3 d_4 s_5  0 ....   0     |
//          |     .       .   .   .  .   ...   0     |
//          |     .          .   .   .               |
//          |     .              .   .   .     0     |
//          |-    0   .   .  0  e_n-2  l_n-1  d_n   -| 
//
//##############################################################################
template<typename T>
class Lower_Triangular_Matrix
{
    std::vector<T> _deltas;    // d
    std::vector<T> _lambdas;   // l
    std::vector<T> _epsilons;  // e
    std::vector<T> _s; 
    T _ZERO = 0.0;
public: 
    Lower_Triangular_Matrix(const int N_hint = 1)
    {
        assert(N_hint > 0);
        _deltas.reserve(N_hint); 
        _lambdas.reserve(N_hint);
        _epsilons.reserve(N_hint);
        _s.reserve(N_hint);
    }
    inline int AddColumn(const T s, const T d, const T l, const T e)
    {
        _s.push_back(s); 
        _deltas.push_back(d); 
        _lambdas.push_back(l); 
        _epsilons.push_back(e); 
        assert(_deltas.size() == _lambdas.size() 
            && _lambdas.size() == _epsilons.size()
            && _s.size() == _deltas.size()); 
        return _deltas.size(); 
    }
    inline T &operator ()(const int row, const int col)
    {
        assert(row>=0 && col>=0);
        const int d = row - col; 
        if      (d== 0) return _deltas.at(col); 
        else if (d== 1) return _lambdas.at(col); 
        else if (d== 2) return _epsilons.at(col);
        else if (d==-1) return _s.at(col); 
        else {_ZERO=(T)0.0; return _ZERO;}
    }
    inline T operator ()(const int row, const int col) const
    {
        const T val = this->operator()(row,col); 
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
void Lower_Triangular_Matrix<T>::
Print()
{
    std::copy(_deltas.begin(), _deltas.end(),
              std::ostream_iterator<T>(std::cout," ")); 
    std::cout << std::endl; 
    std::copy(_lambdas.begin(), _lambdas.end(),
              std::ostream_iterator<T>(std::cout," ")); 
    std::cout << std::endl; 
    std::copy(_epsilons.begin(), _epsilons.end(),
              std::ostream_iterator<T>(std::cout," ")); 
    std::cout << std::endl; 
}

template<typename T>
void Lower_Triangular_Matrix<T>::
Print(const int N)
{
    Print(N,N);
}

template<typename T>
void Lower_Triangular_Matrix<T>::
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
