#ifndef SPARSE_MATRIX_HPP
#define SPARSE_MATRIX_HPP
#include "Eigen/Dense"
//##############################################################################
// Class Sparse_Matrix
//   Custom sparse matrix wrapper. Need to support 
//   A v and A^* v matrix-vector operations.
//##############################################################################
template<typename T, int N> 
class Sparse_Matrix
{
public: 
    using USYM_Vector = Eigen::Matrix<T,N,1>; 
    using USYM_VectorX= Eigen::Matrix<T,Eigen::Dynamic,1>; 
    using USYM_Matrix = Eigen::Matrix<T,N,N>; 
    using USYM_MatrixX= Eigen::Matrix<T,N,Eigen::Dynamic>; 
private: 
    USYM_Matrix _m; 
public:
    Sparse_Matrix() = default; 
    Sparse_Matrix(const USYM_Matrix &m)
        : _m(std::move(m))
    {}
    Sparse_Matrix(Sparse_Matrix &&rhs)
        : _m(std::move(rhs._m))
    {}
    Sparse_Matrix &operator=(Sparse_Matrix &&rhs)
    {
        _m = std::move(rhs._m);
    }
    USYM_Vector operator*(const USYM_Vector &rhs) const
    {
        return _m*rhs;
    }
    void Premultiply_By_Matrix(const USYM_Vector &v_i, USYM_Vector &v_o); 
    void Premultiply_By_Matrix_Conjugate(const USYM_Vector &v_i, 
                                         USYM_Vector &v_o); 
};

//############################################################################## 
// Function Premultiply_By_Matrix
//##############################################################################
template<typename T, int N> 
void Sparse_Matrix<T,N>::
Premultiply_By_Matrix(const USYM_Vector &v_i, USYM_Vector &v_o)
{
    v_o = _m * v_i; 
}

//##############################################################################
// Function Premultiply_By_Matrix_Conjugate
//##############################################################################
template<typename T, int N> 
void Sparse_Matrix<T,N>::
Premultiply_By_Matrix_Conjugate(const USYM_Vector &v_i, USYM_Vector &v_o)
{
    v_o = _m.transpose() * v_i; 
}

#endif
