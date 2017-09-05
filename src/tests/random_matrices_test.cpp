#include <memory>
#include <fstream>
#include <iostream>
#include "macros.h"
#include "usym_linear_solver.hpp"
#include "sparse_matrix.hpp"
#include "usym_tridiag.hpp"
#include "Eigen/Dense"
//#define GENERATE_COMPATIBLE_SYSTEM
//##############################################################################
using T = double; 
using T_Vector = Eigen::Matrix<T,Eigen::Dynamic,1>; 
using T_Matrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>; 
//##############################################################################
inline bool File_Exists(const std::string &name) {
    std::ifstream f(name.c_str());
    return f.good();
}
//##############################################################################
int main(int argc, char **argv)
{
    if (argc != 5) 
    {
        std::cerr << "**Usage: " << argv[0] 
                  << " <mode: 0 for USYMLQ; 1 for USYMQR> <M> <N> <data_dir>\n";
        exit(1); 
    }
    const auto mode = (atof(argv[1]) == 0 ? USYMLQ : USYMQR); 
    const int M = atoi(argv[2]);
    const int N = atoi(argv[3]);
    const std::string datadir(argv[4]);

    const std::string modeName = (mode == USYMLQ ? "USYMLQ" : "USYMQR");
    const int maxStep = std::max(M,N)*5;
    // initialize A and b
    std::shared_ptr<T_Matrix> A(new T_Matrix);
    (*A) = T_Matrix::Random(M,N); 
#ifdef GENERATE_COMPATIBLE_SYSTEM
    const T_Vector xstar = T_Vector::Random(N); 
    const T_Vector b  = (*A)*xstar; 
#else
    const T_Vector b = T_Vector::Random(M);
#endif
    // output storage
    T_Vector x; 
    T rnorm; 

    // determine filename
    std::string filename_out, filename_A, filename_b, filename_x;
    filename_out = datadir + "/" + modeName + ".log"; 
    filename_A   = datadir + "/" + "A.txt"; 
    filename_b   = datadir + "/" + "b.txt"; 
    filename_x   = datadir + "/" + modeName + "_x.txt";
    if (File_Exists(filename_out))
    {
        std::cout << "**WARNING** Log file exists, abort solve: " 
                  << filename_out << std::endl; 
        exit(1);
    }
    std::ofstream file(filename_out.c_str()); 
    // initialize solver
    USYM_Linear_Solver<T,T_Vector,T_Matrix> solver(A,b); 
    solver.Initialize(); 
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
    file.open(filename_x.c_str(), std::ios::out); 
    file.precision(16);
    file << x << std::endl;
    file.close(); 
}
