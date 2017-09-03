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
inline bool File_Exists(const char *name) {
    std::ifstream f(name);
    return f.good();
}
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

    const std::string modeName = (mode == USYMLQ ? "USYMLQ" : "USYMQR");
    const int maxStep = std::max(M,N)*100;
    // initialize A and b
    std::shared_ptr<T_Matrix> A(new T_Matrix);
    (*A) = T_Matrix::Random(M,N); 
    const T_Vector b = T_Vector::Random(M); 
    const T_Vector x0 = T_Vector::Zero(N); 
    // output storage
    T_Vector x; 
    T rnorm; 

    // determine filename
    std::string filename_out, filename_A, filename_b;
    { 
        char buf[512]; 
        int count = 0;
        while(true)
        {
            snprintf(buf, 512, 
                     "data/random_mat/log_%d_M%d_N%d_%s.log",
                     count, M, N, modeName.c_str());
            if (!File_Exists(&(buf[0])))
            {
                filename_out = std::string(buf);
                snprintf(buf, 512, 
                         "data/random_mat/A_%d_M%d_N%d_%s.txt",
                         count, M, N, modeName.c_str());
                filename_A   = std::string(buf);
                snprintf(buf, 512, 
                         "data/random_mat/b_%d_M%d_N%d_%s.txt",
                         count, M, N, modeName.c_str());
                filename_b   = std::string(buf);
                break;
            }
            count ++;
        }
    }
    std::ofstream file(filename_out.c_str()); 
    // initialize solver with x0
    USYM_Linear_Solver<T,T_Vector,T_Matrix> solver(A,b); 
    solver.Initialize(x0); 
    solver.Set_Mode(mode);
    solver.SetMaxIteration(maxStep);
    solver.Set_Logging(&(file)); 
    solver.Set_Verbose_Level(2);
    solver.Solve(x, rnorm); 
    file.close();

    file.open(filename_A.c_str(), std::ios::out); 
    file.precision(16);
    file << (*A) << std::endl;
    file.close(); 
    file.open(filename_b.c_str(), std::ios::out); 
    file.precision(16);
    file << b << std::endl;
    file.close(); 
}
