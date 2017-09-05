import numpy as np
import pickle
################################################################################
## Class Solver_Results
################################################################################
class Solver_Results: 
    def __init__(self): 
        self.stype    = None
        self.maxItn   = None
        self.a_tol    = None
        self.b_tol    = None
        self.verbose  = None
        self.mat_size = None

        self.steps    = None
        self.itnv     = None
        self.x1       = None
        self.resv     = None
        self.normA    = None
        self.relNormAr= None
        self.relRes   = None

        self.timing   = None

        self.flag     = None
        self.r_exact  = None

    @staticmethod
    def PickleFile(prefix): 
        return prefix + '.dat'

    @staticmethod
    def Save(results, save_file): 
        with open(save_file, 'wb') as stream:
            pickle.dump(results, stream, 2)

    @staticmethod
    def Load(load_file): 
        data = None
        with open(load_file, 'rb') as stream:
            data = pickle.load(stream)
        return data

    @staticmethod
    def Parse(filename): 
        lines = open(filename,'r').readlines()
        results = None
        # parse header
        ii = 0
        while ii < (len(lines)): 
            l = lines[ii]
            if l.find('Solver Header') != -1: 
                results = Solver_Results()
                while True:
                    ii += 1
                    l = lines[ii]
                    if l.find('==========') != -1: 
                        break
                    elif l.find('solver type'  ) != -1: 
                        results.stype = l.split(':')[-1].split(' ')[-1][:-1]
                    elif l.find('max iteration') != -1: 
                        results.maxItn = int(l.split(':')[-1].split(' ')[-1])
                    elif l.find('a_tol'        ) != -1: 
                        results.a_tol = float(l.split(':')[-1].split(' ')[-1])
                    elif l.find('b_tol'        ) != -1: 
                        results.b_tol = float(l.split(':')[-1].split(' ')[-1])
                    elif l.find('verbose level') != -1: 
                        results.verbose = int(l.split(':')[-1].split(' ')[-1])
                    elif l.find('matrix size'  ) != -1: 
                        s = l.split(':')[-1].split(' ')[-1]
                        s = s.split('-by-')
                        results.mat_size = [int(s[0]), int(s[1])]
            elif l.find('Solver START') != -1: 
                steps = 0
                jj = ii
                jj += 2
                while True:
                    jj += 1
                    l = lines[jj]
                    if l.find('Solver END') != -1: 
                        break
                    steps += 1
                results.steps    = steps
                results.itnv     = np.zeros(steps, dtype=int  )
                results.x1       = np.zeros(steps, dtype=float)
                results.resv     = np.zeros(steps, dtype=float)
                results.normA    = np.zeros(steps, dtype=float)
                results.relNormAr= np.zeros(steps, dtype=float)
                results.relRes   = np.zeros(steps, dtype=float)
                ii += 2
                steps = 0
                while True: 
                    ii += 1
                    l = lines[ii]
                    tokens = l.split()
                    results.itnv     [steps] =   int(tokens[0])
                    results.x1       [steps] = float(tokens[1])
                    results.resv     [steps] = float(tokens[2])
                    results.normA    [steps] = float(tokens[3])
                    results.relNormAr[steps] = float(tokens[4])
                    results.relRes   [steps] = float(tokens[5])
                    steps += 1
                    if steps >= results.steps-1:
                        break
            elif l.find('Timing statistics') != -1: 
                ii += 1
                l = lines[ii]
                results.timing = l.split(':')[-1].split(' ')[1]
            elif l.find('Solver Footer') != -1: 
                ii += 1
                l = lines[ii]
                results.flag    =   int(l.split(':')[-1].split(' ')[1])
                ii += 3
                l = lines[ii]
                results.r_exact = float(l.split(':')[-1].split(' ')[1])
            ii += 1
        return results
