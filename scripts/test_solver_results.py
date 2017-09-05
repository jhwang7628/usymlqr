#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt 

# stop flags
flags0 = sorted(np.loadtxt('data/stopflags0.txt', dtype=int))
flags1 = sorted(np.loadtxt('data/stopflags1.txt', dtype=int))

flags_count0 = dict()
flags_count1 = dict()
for ii in flags0: 
    if ii in flags_count0: 
        flags_count0[ii] += 1
    else: 
        flags_count0[ii] = 0

for ii in flags1: 
    if ii in flags_count1: 
        flags_count1[ii] += 1
    else: 
        flags_count1[ii] = 0

print 'Flag count for USYMLQ:', flags_count0
print 'Flag count for USYMQR:', flags_count1


# final residual
resfinal0 = sorted(np.loadtxt('data/resfinal0.txt', dtype=float))
resfinal1 = sorted(np.loadtxt('data/resfinal1.txt', dtype=float))
bin_edges   = np.logspace(-17, 0, 18)
bin_centers = np.linspace(-16.5, 0.5, 17)
print bin_centers
hist0 = np.histogram(resfinal0, bin_edges)
hist1 = np.histogram(resfinal1, bin_edges)
print hist0
print hist1
plt.figure()
plt.bar(hist0[0])
plt.figure()
plt.bar(hist1[0])
plt.show()
