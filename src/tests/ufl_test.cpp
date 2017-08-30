#include <memory>
#include <fstream>
#include <iostream>
#include "macros.h"
#include "usym_linear_solver.hpp"
#include "sparse_matrix.hpp"
#include "usym_tridiag.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "ufl_io/mmio.h"
//##############################################################################
using T = double; 
using T_Vector = Eigen::Matrix<T,Eigen::Dynamic,1>; 
using T_Matrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>; 
//##############################################################################
int main(int argc, char **argv)
{
    if (argc != 3) 
    {
        std::cerr << "**Usage: " << argv[0] 
                  << " <mode: 0 for USYMLQ; 1 for USYMQR> <matrix-market-filename>\n";
        exit(1); 
    }
    const auto mode = (atof(argv[1]) == 0 ? USYMLQ : USYMQR); 
    const char *filename = argv[2]; 

    // error checking
    FILE *file; 
    MM_typecode matcode; 
    if ((file = fopen(filename, "r")) == nullptr)
    {
        std::cerr << "**ERROR** Cannot open file: " << filename << std::endl;
        exit(2); 
    }
    if (mm_read_banner(file, &matcode) != 0) 
    {
        std::cerr << "**ERROR** Could not process Matrix Market banner.\n";
        exit(3); 
    }

    // initialize A and b
    int M, N, nz, retcode; 
    if ((retcode = mm_read_mtx_crd_size(file, &M, &N, &nz)) != 0) 
        exit(4); 
    std::cout << M << " " << N << " " << nz << std::endl;

#if 0
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
#endif
}
