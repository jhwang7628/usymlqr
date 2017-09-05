# USYMLQ/USYMQR
USYMLQ and USYMQR algorithms generalize SYMMLQ and MINRES for solving
large, sparse, unsymmetric linear system of equations. 
In the symmetric case, USYMLQ falls back to SYMMLQ and USYMQR falls 
back to MINRES. It is based on an efficient orthogonal tridiagonalization 
procedure described in Saunders et al. 1988 paper 
"Two Conjugate-Gradient-Type Methods For Unsymmetric Linear Equations". 
USYMLQ can be used for square or under-determined systems, whereas
USYMQR can be used for square, under-determined, and over-determined
systems. They are iterative solvers and each step has O(N) operations;
the overall storage is O(N) as well (unlike GMRES). 

The solver was implemented in templated C++ with Eigen as backend 
linear algebra library. It has no external dependency otherwise. 

This is part of my final project for CME338 at Stanford, taught 
by Michael Saunders. 

## Build Instructions
To build the project, simply run
```
make
```

## Usage
The solver encapsulates the tridiagonalization and linear solve. 
Below is an example of how to use it:
```
// ... define A and b and ...
T_Vector x;
T rnorm;
USYM_Linear_Solver<T,T_Vector,T_Matrix> solver(A,b);
solver.Initialize();
solver.Set_Mode(USYMQR);
solver.SetMaxIteration(2000);
solver.Set_Tol(1E-12, 1E-12);
solver.Solve(x, rnorm);
```

## Documentation
Please see doc/documentation.pdf for more detailed analysis and 
solver testing.
