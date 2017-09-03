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
    M = [4, 10, 50, 100, 500, 1000, 2500, 5000]
    N = [4, 10, 50, 100, 500, 1000, 2500, 5000]
    # for ii in range(len(M)): 
    #     m = M[ii]
    #     n = N[ii]
    #     cmd = 'bin/random_matrices_test 0 %d %d' %(m,n)
    #     Call_Cmd(cmd)
    #     cmd = 'bin/random_matrices_test 1 %d %d' %(m,n)
    #     Call_Cmd(cmd)

    M = M
    N = [x/2 for x in N]
    for ii in range(len(M)): 
        m = M[ii]
        n = N[ii]
        cmd = 'bin/random_matrices_test 0 %d %d' %(m,n)
        Call_Cmd(cmd)
        cmd = 'bin/random_matrices_test 1 %d %d' %(m,n)
        Call_Cmd(cmd)

if __name__ == "__main__": 
    Run_Random()

