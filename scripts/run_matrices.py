#!/usr/bin/env python 
import glob
from subprocess import call

collection = 'matrices/list_maxrows_10000_maxcols_10000_real'
# collection = 'matrices/tests'
groups = glob.glob('%s/*' %(collection))
for g in groups: 
    datadir = 'data/%s' %(g)
    call('mkdir -p %s' %(datadir), shell=True)
    mtxs = glob.glob('%s/*.mtx' %(g))
    for mtxfile in mtxs: 
        basename = mtxfile.split('/')[-1]
        logfile = '%s/%s.log' %(datadir, basename[:-4])
        cmd = 'bin/ufl_test %s 2>&1 > %s' %(mtxfile, logfile)
        print cmd
        call(cmd, shell=True)


