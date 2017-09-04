#!/usr/bin/env python 
import sys,os
import numpy as np
import matplotlib.pyplot as plt
from solver_results        import *
from matlab_solver_results import *

data_dir = 'data/random_mat/benchmark_compatible_5000_2500'
if data_dir.split('/')[-1].split('_') == 2: square = True # hack..
else                                      : square = False

if not os.path.isdir(data_dir):
    sys.exit()

A_file      = '%s/A.txt'      %(data_dir)
b_file      = '%s/b.txt'      %(data_dir)
USYMLQ_file = '%s/USYMLQ.log' %(data_dir)
USYMQR_file = '%s/USYMQR.log' %(data_dir)

A_svs_file  = '%s/A_singular_values.txt' %(data_dir)
A_svs       = np.loadtxt(A_svs_file)
if not square: 
    lm = 30

print 'Reading'
usymlq = Solver_Results.Parse(USYMLQ_file)
usymqr = Solver_Results.Parse(USYMQR_file)
lsqrma = Matlab_Solver_Results.Parse(data_dir, 'lsqr')
lsmrma = Matlab_Solver_Results.Parse(data_dir, 'lsmr')
if square: 
    bcgsma = Matlab_Solver_Results.Parse(data_dir, 'bcgs')
    bcgspm = Matlab_Solver_Results.Parse(data_dir, 'bcgsp')

print 'Plotting'
plt.figure(figsize=[9,6])
plt.semilogy(usymlq.itnv, usymlq.resv, '-', linewidth=1.5, label='USYMLQ')
plt.semilogy(usymqr.itnv, usymqr.resv, '-', linewidth=1.5, label='USYMQR')
plt.semilogy(lsqrma.itnv, lsqrma.resv, '-', linewidth=1.5, label='LSQR')
plt.semilogy(lsmrma.itnv, lsmrma.resv, '-', linewidth=1.5, label='LSMR')
if square:
    plt.semilogy(bcgsma.itnv, bcgsma.resv, linewidth=1.5, label='BiCGStab')
    plt.semilogy(bcgspm.itnv, bcgspm.resv, linewidth=1.5, label='BiCGStab-ILU')
    if bcgsma.flag == 4: # unsuccessful
        plt.xlim([0, max(max(len(usymlq.resv),
                             len(usymqr.resv)),len(lsqrma.resv))])
plt.xlabel('Iteration')
plt.ylabel('norm(b - Ax_i)')
plt.legend(loc=0)
plt.savefig('%s/convergence.pdf' %(data_dir))

plt.figure(figsize=[9,6])
plt.semilogy(usymlq.itnv[:lm], usymlq.resv[:lm], 'o-', linewidth=1.5, label='USYMLQ')
plt.semilogy(usymqr.itnv[:lm], usymqr.resv[:lm], '^-', linewidth=1.5, label='USYMQR')
plt.semilogy(lsqrma.itnv[:lm], lsqrma.resv[:lm], 's-', linewidth=1.5, label='LSQR')
plt.semilogy(lsmrma.itnv[:lm], lsmrma.resv[:lm], '*-', linewidth=1.5, label='LSMR')
if square:
    plt.semilogy(bcgsma.itnv[:lm], bcgsma.resv[:lm], 'x-', linewidth=1.5, label='BiCGStab')
    plt.semilogy(bcgspm.itnv[:lm], bcgspm.resv[:lm], '>-', linewidth=1.5, label='BiCGStab-ILU')
    if bcgsma.flag == 4: # unsuccessful
        plt.xlim([0, max(max(len(usymlq.resv),
                             len(usymqr.resv)),len(lsqrma.resv))])
plt.xlabel('Iteration')
plt.ylabel('norm(b - Ax_i)')
plt.legend(loc=0)
plt.savefig('%s/convergence_init.pdf' %(data_dir))

plt.figure(figsize=[6,6])
plt.plot(A_svs[::-1], '-o')
plt.title('Singular Values of A')
plt.xlabel('Number')
plt.ylabel('Singular Values')
plt.savefig('%s/A_singular_values.pdf' %(data_dir))
plt.show()
