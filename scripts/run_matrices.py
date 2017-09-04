#!/usr/bin/env python 
import glob,os
import numpy as np
from subprocess import call

def Call_Cmd(cmd): 
    print cmd
    call(cmd, shell=True)

def Run_UFL():
    collection = 'matrices/list_maxrows_10000_maxcols_10000_real'
    # collection = 'matrices/tests'
    groups = glob.glob('%s/*' %(collection))
    for g in groups: 
        datadir = 'data/%s' %(g)
        call('mkdir -p %s' %(datadir), shell=True)
        mtxs = glob.glob('%s/*.mtx' %(g))
        for mtxfile in mtxs: 
            basename = mtxfile.split('/')[-1]
            logfile = '%s/%s_USYMQR.log' %(datadir, basename[:-4])
            if os.path.isfile(logfile): 
                continue
            cmd = 'bin/ufl_test %s 1 2>&1 > %s' %(mtxfile, logfile)
            Call_Cmd(cmd)

def Run_Random(): 
    M = [50, 500, 5000]
    N = [x/2 for x in M]
    for ii in range(len(M)): 
        m = M[ii]
        n = N[ii]
        datadir = 'data/random_mat/benchmark_%d_%d' %(m,n)
        if not os.path.isdir(datadir): 
            Call_Cmd('mkdir -p %s' %(datadir))
        Call_Cmd('bin/random_matrices_test 0 %d %d %s' %(m,n, datadir))
        Call_Cmd('bin/random_matrices_test 1 %d %d %s' %(m,n, datadir))

def Run_Compatible_Random(): 
    M = [50, 500, 5000]
    N = [x/2 for x in M]
    for ii in range(len(M)): 
        m = M[ii]
        n = N[ii]
        datadir = 'data/random_mat/benchmark_compatible_%d_%d' %(m,n)
        if not os.path.isdir(datadir): 
            Call_Cmd('mkdir -p %s' %(datadir))
        Call_Cmd('bin/random_matrices_test 0 %d %d %s' %(m,n, datadir))
        Call_Cmd('bin/random_matrices_test 1 %d %d %s' %(m,n, datadir))

if __name__ == "__main__": 
    # Run_Random()
    Run_Compatible_Random()

