#!/usr/bin/env python 
import sys,os
import numpy as np
import matplotlib.pyplot as plt
from solver_results        import *
from matlab_solver_results import *

data_dir = 'data/random_mat/benchmark_50_25'
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
usymqr = Solver_Results.Parse(USYMQR_file)
lsqrma = Matlab_Solver_Results.Parse(data_dir, 'lsqr')
lsmrma = Matlab_Solver_Results.Parse(data_dir, 'lsmr')
if square: 
    usymlq = Solver_Results.Parse(USYMLQ_file)
    bcgsma = Matlab_Solver_Results.Parse(data_dir, 'bcgs')
    bcgspm = Matlab_Solver_Results.Parse(data_dir, 'bcgsp')

print 'Plotting'
plt.figure(figsize=[9,6])
if square:
    plt.semilogy(usymlq.itnv, usymlq.resv, '-', linewidth=1.5, label='USYMLQ')
    plt.semilogy(usymqr.itnv, usymqr.resv, '-', linewidth=1.5, label='USYMQR')
    plt.semilogy(lsqrma.itnv, lsqrma.resv, '-', linewidth=1.5, label='LSQR')
    plt.semilogy(lsmrma.itnv, lsmrma.resv, '-', linewidth=1.5, label='LSMR')
    plt.semilogy(bcgsma.itnv, bcgsma.resv, linewidth=1.5, label='BiCGStab')
    plt.semilogy(bcgspm.itnv, bcgspm.resv, linewidth=1.5, label='BiCGStab-ILU')
    if bcgsma.flag == 4: # unsuccessful
        plt.xlim([0, max(max(len(usymlq.resv),
                             len(usymqr.resv)),len(lsqrma.resv))])
    plt.ylabel('norm(b - Ax_i)')
else: 
    plt.semilogy(usymqr.itnv[1:], usymqr.relNormAr[1:], '-', linewidth=1.5, label='USYMQR')
    plt.semilogy(lsqrma.itnv[1:], lsqrma.relNormAr[1:], '-', linewidth=1.5, label='LSQR')
    plt.semilogy(lsmrma.itnv[1:], lsmrma.relNormAr[1:], '-', linewidth=1.5, label='LSMR')
    plt.ylabel('norm(A\'(b - Ax_i))/norm(A)')
plt.xlabel('Iteration')
plt.legend(loc=0)
plt.savefig('%s/convergence.pdf' %(data_dir))

plt.figure(figsize=[6,6])
plt.plot(A_svs[::-1], '-o')
plt.title('Singular Values of A')
plt.xlabel('Number')
plt.ylabel('Singular Values')
plt.savefig('%s/A_singular_values.pdf' %(data_dir))
plt.show()
