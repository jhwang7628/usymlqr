#!/usr/bin/env python 
import glob,os
import numpy as np
from subprocess import call
from Queue import Queue
from threading import Thread, current_thread

def Call_Cmd(cmd): 
    print cmd
    call(cmd, shell=True)

def job(q,tid):
    while True: 
        cmd = q.get()
        Call_Cmd(cmd)
        q.task_done()

def Run_UFL():
    mode = 'USYMQR'
    num_threads = 48

    if mode == 'USYMQR': 
        modeName = 'USYMQR'
        modeId   = 1
    elif mode == 'USYMLQ': 
        modeName = 'USYMLQ'
        modeId   = 0
    else: 
        return

    collection = 'matrices/list_maxrows_10000_maxcols_10000_real'
    q = Queue(maxsize=0)

    for tid in range(num_threads):
        worker = Thread(target=job, args=(q,tid,))
        worker.setDaemon(True)
        worker.start()

    # collection = 'matrices/tests'
    groups = glob.glob('%s/*' %(collection))
    for g in groups: 
        datadir = 'data/%s' %(g)
        call('mkdir -p %s' %(datadir), shell=True)
        mtxs = glob.glob('%s/*.mtx' %(g))
        for mtxfile in mtxs: 
            basename = mtxfile.split('/')[-1]
            logfile = '%s/%s_%s.log' %(datadir, basename[:-4], modeName)
            if os.path.isfile(logfile): 
                continue
            cmd = 'bin/ufl_test %s %d 2>&1 > %s' %(mtxfile, modeId, logfile)
            q.put(cmd)
            # Call_Cmd(cmd)

    q.join()

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
        datadir = 'data/random_mat/benchmark_%d_%d' %(m,n)
        if not os.path.isdir(datadir): 
            Call_Cmd('mkdir -p %s' %(datadir))
        # Call_Cmd('bin/random_matrices_test 0 %d %d %s' %(m,n, datadir))
        Call_Cmd('bin/random_matrices_test 1 %d %d %s' %(m,n, datadir))

if __name__ == "__main__": 
    # Run_Random()
    # Run_Compatible_Random()
    Run_UFL()

