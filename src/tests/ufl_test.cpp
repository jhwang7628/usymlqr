#include <memory>
#include <fstream>
#include <iostream>
#include "macros.h"
#include "usym_linear_solver.hpp"
#include "sparse_matrix.hpp"
#include "usym_tridiag.hpp"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "ufl_io/mmio.h"
//##############################################################################
using T = double; 
using T_Vector = Eigen::Matrix<T,Eigen::Dynamic,1>; 
using T_Matrix = Eigen::SparseMatrix<T>; 
//##############################################################################
int main(int argc, char **argv)
{
    if (argc < 2) 
    {
        std::cout << "**Usage: " << argv[0] 
                  << " <matrix-market-filename> [mode: 0 for USYMLQ; 1 for USYMQR; default: USYMQR]\n";
        exit(1); 
    }
    const char *filename = argv[1]; 

    // error checking
    FILE *file; 
    MM_typecode matcode; 
    if ((file = fopen(filename, "r")) == nullptr)
    {
        std::cout << "**ERROR** Cannot open file: " << filename << std::endl;
        exit(2); 
    }
    if (mm_read_banner(file, &matcode) != 0) 
    {
        std::cout << "**ERROR** Could not process Matrix Market banner.\n";
        exit(3); 
    }
    if ((!mm_is_matrix    (matcode)) || 
        (!mm_is_coordinate(matcode)) || 
        (!mm_is_real      (matcode)) || 
        (!mm_is_sparse    (matcode)) || 
        (!mm_is_general   (matcode)))
    {
        std::cout << "**ERROR** Matrix Market format not supported: "
                  << mm_typecode_to_str(matcode) << std::endl;
        exit(4);
    }

    // initialize A and b
    int M, N, nz, retcode; 
    if ((retcode = mm_read_mtx_crd_size(file, &M, &N, &nz)) != 0) 
        exit(5); 

    // default to QR if not specified
    auto mode = (argc == 3 ? (atof(argv[2]) == 0 ? USYMLQ : USYMQR) : USYMQR); 
        
    std::cout << "Reading matrix: " << M << "-by-" << N << " with " 
              << nz << " non-zero entries\n"; 
    std::shared_ptr<T_Matrix> A(new T_Matrix(M,N));
    int buf_i, buf_j; 
    double buf_v; 
    for (int ii=0; ii<nz; ++ii)
    {
        fscanf(file, "%d %d %lg\n", &buf_i, &buf_j, &buf_v); 
        -- buf_i; // convert to 0-based 
        -- buf_j;
        A->insert(buf_i, buf_j) = (T)buf_v; 
    }
    fclose(file); 

    // settings
    const int maxStep = std::max(M,N)*5;
    // constructing b
    T_Vector x0    = T_Vector::Zero(N); 
    T_Vector xstar = T_Vector::Random(N); 
    T_Vector b     = (*A)*xstar; 
    // solution vector
    T_Vector x; 
    T rnorm; 

    // initialize solver with x0
    USYM_Linear_Solver<T,T_Vector,T_Matrix> solver(A,b); 
    solver.Initialize(x0); 
    solver.Set_Exact_Solution(xstar); 
    solver.Set_Mode(mode);
    solver.SetMaxIteration(maxStep);
    solver.Set_Verbose_Level(2);
    solver.Solve(x, rnorm); 
}
