import glob
import numpy as np
################################################################################
## Class Matlab_Solver_Results
################################################################################
class Matlab_Solver_Results: 
    def __init__(self): 
        self.data_dir = None
        self.stype    = None

        self.steps    = None
        self.resv     = None
        self.timing   = None
        self.flag     = None
        self.r_exact  = None

        self.itnv     = None

        # optional
        self.maxItn   = None
        self.a_tol    = None
        self.mat_size = None

    @staticmethod
    def Parse(data_dir, solver_type): 
        # solver_type: lsqr | lsmr | bcgs | bcgsp 
        results = Matlab_Solver_Results() 
        results.data_dir = data_dir
        results.stype    = solver_type
        d = data_dir; s = solver_type

        results.steps  = int(open('%s/%s_ite.txt'  %(d, s)).readline(1))
        results.flag   = int(open('%s/%s_flag.txt' %(d, s)).readline(1))
        results.resv   = np.loadtxt('%s/%s_resv.txt' %(d, s))
        results.timing = np.loadtxt('%s/%s_t.txt'    %(d, s))
        results.r_exact= np.loadtxt('%s/%s_res.txt'  %(d, s))

        results.itnv   = range(len(results.resv))

        return results
