#!/usr/bin/env python
from solver_results import *

datafile = '/Users/jui-hsien/code/usymqr/data/random_mat/benchmark_50/USYMLQ.log'
testfile = 'tmp.pkl'
results  = Solver_Results.Parse(datafile)
Solver_Results.Save(results, testfile)
results2 = Solver_Results.Load(testfile)

assert(results.stype == results2.stype)
assert(results.maxItn == results2.maxItn)
assert(results.a_tol == results2.a_tol)
assert(results.b_tol == results2.b_tol)
assert(results.verbose == results2.verbose)
assert(results.mat_size == results2.mat_size)
assert(results.steps == results2.steps)
assert(results.itnv.all() == results2.itnv.all())
assert(results.x1.all() == results2.x1.all())
assert(results.resv.all() == results2.resv.all())
assert(results.normA.all() == results2.normA.all())
assert(results.relNormAr.all() == results2.relNormAr.all())
assert(results.relRes.all() == results2.relRes.all())
assert(results.timing == results2.timing)
assert(results.flag == results2.flag)
assert(results.r_exact == results2.r_exact)

print 'Passed'
