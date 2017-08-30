#include <memory>
#include <fstream>
#include <iostream>
#include "macros.h"
#include "usym_linear_solver.hpp"
#include "sparse_matrix.hpp"
#include "usym_tridiag.hpp"
#include "Eigen/Dense"
//##############################################################################
using T = double; 
using T_Vector = Eigen::Matrix<T,Eigen::Dynamic,1>; 
using T_Matrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>; 
//##############################################################################
int main(int argc, char **argv)
{
    if (argc != 4) 
    {
        std::cerr << "**Usage: " << argv[0] 
                  << " <mode: 0 for USYMLQ; 1 for USYMQR> <M> <N>\n";
        exit(1); 
    }
    const auto mode = (atof(argv[1]) == 0 ? USYMLQ : USYMQR); 
    const int M = atoi(argv[2]);
    const int N = atoi(argv[3]);

    const int maxStep = std::max(M,N)*100;
    // initialize A and b
    std::shared_ptr<T_Matrix> A(new T_Matrix);
    (*A) = T_Matrix::Random(M,N); 
    const T_Vector b = T_Vector::Random(M); 
    const T_Vector x0 = T_Vector::Zero(N); 
    // output storage
    T_Vector x; 
    T rnorm; 
    std::ofstream file("tmp.txt"); 
    // initialize solver with x0
    USYM_Linear_Solver<T,T_Vector,T_Matrix> solver(A,b); 
    solver.Initialize(x0); 
    solver.Set_Mode(mode);
    solver.SetMaxIteration(maxStep);
    solver.Set_Logging(&(std::cout)); 
    solver.Set_Verbose_Level(2);
    solver.Solve(x, rnorm); 
    file.close();

    file.open("A.txt", std::ios::out); 
    file << (*A) << std::endl;
    file.close(); 
    file.open("b.txt", std::ios::out); 
    file << b << std::endl;
    file.close(); 
}
