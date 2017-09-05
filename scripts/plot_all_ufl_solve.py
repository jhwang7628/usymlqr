#!/usr/bin/env python
import sys,os
import numpy as np 
import matplotlib.pyplot as plt
from solver_results import *

##
listfile = sys.argv[1]
##
mat_sizes0 = []
mat_sizes1 = []
timings0   = []
timings1   = []
stopflags0 = []
stopflags1 = []
resreduce0 = []
resreduce1 = []
with open(listfile,'r') as stream: 
    logs = stream.readlines()
    for logfile in logs: 
        logfile = logfile[:-1]
        # print 'reading %s' %(logfile.split('/')[-1])
        if not os.path.isfile(logfile): 
            continue

        # determine if should read cache or read from log
        cache_file = Solver_Results.PickleFile(logfile[:-4])
        cache = False
        if os.path.isfile(cache_file): 
            results = Solver_Results.Load(cache_file)
        else: 
            try:
                results = Solver_Results.Parse(logfile)
                cache = True
            except IndexError:  # file not intact
                continue
        if results is None: 
            continue
        if cache: Solver_Results.Save(results, cache_file)

        M,N = results.mat_size
        if N > 10000: 
            continue

        s = results.mat_size
        if   results.stype == 'USYMLQ': 
            mat_sizes0.append(s)
            timings0.append(results.timing)
            stopflags0.append(results.flag)
            resreduce0.append(results.resv[-2]/results.resv[0])
        elif results.stype == 'USYMQR': 
            mat_sizes1.append(s)
            timings1.append(results.timing)
            stopflags1.append(results.flag)
            resreduce1.append(results.resv[-2]/results.resv[0])

        if (len(resreduce1)>0 and resreduce1[-1] > 1E-2) or (len(resreduce0)>0 and resreduce0[-1] > 1E-2): 
            print M,N,logfile
        
    mat_sizes0 = np.array(mat_sizes0, dtype=int)
    mat_sizes1 = np.array(mat_sizes1, dtype=int)
    
    # matrix size
    plt.figure(figsize=[6,6])
    print mat_sizes1
    plt.plot(mat_sizes1[:,0], mat_sizes1[:,1], 'ob', label='USYMQR')
    plt.plot(mat_sizes0[:,0], mat_sizes0[:,1], 'xr', label='USYMLQ')
    plt.xlabel('M')
    plt.ylabel('N')
    plt.axis('equal')
    plt.legend(loc=4)
    plt.savefig('doc/figures/ufl_problems_10k.pdf')

    # termination time
    MN0 = np.multiply(mat_sizes0[:,0], mat_sizes0[:,1])
    MN1 = np.multiply(mat_sizes1[:,0], mat_sizes1[:,1])
    plt.figure(figsize=[9,6])
    plt.loglog(MN1, timings1, 'ob', label='USYMQR')
    plt.loglog(MN0, timings0, 'xr', label='USYMLQ')
    plt.xlabel('Problem size (MN)')
    plt.ylabel('Termination time (sec)')
    plt.legend(loc=2)
    plt.savefig('doc/figures/ufl_timing_10k.pdf')

    # flag
    plt.figure(figsize=[9,6])
    plt.semilogx(MN1, stopflags1, 'ob', label='USYMQR')
    plt.semilogx(MN0, stopflags0, 'xr', label='USYMLQ')
    plt.xlabel('Problem size (MN)')
    plt.ylabel('Termination flag')
    plt.legend(loc=0)
    plt.savefig('doc/figures/ufl_flag_10k.pdf')

    # residual reduction
    plt.figure(figsize=[9,6])
    plt.loglog(MN1, resreduce1, 'ob', label='USYMQR')
    plt.loglog(MN0, resreduce0, 'xr', label='USYMLQ')
    plt.xlabel('Problem size (MN)')
    plt.ylabel('Residual reduction norm(r)/norm(b)')
    plt.legend(loc=2)
    plt.savefig('doc/figures/ufl_residual_reduction_10k.pdf')

    plt.show()

