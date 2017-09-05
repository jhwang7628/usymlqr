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
bin_edges   = np.logspace(-16, 0, 17)
bin_centers = np.linspace(-15.5, -0.5, 16)
hist0 = np.histogram(resfinal0, bin_edges)
hist1 = np.histogram(resfinal1, bin_edges)
plt.figure()
ax = plt.subplot(111)
ax.bar(bin_centers+0.175, hist1[0], width=0.35, color='b', align='center', label='USYMQR')
ax.bar(bin_centers-0.175, hist0[0], width=0.35, color='r', align='center', label='USYMLQ')
ax.set_xticks(np.log10(bin_edges[1:]))
ax.set_xlabel('log(norm(r))')
ax.set_ylabel('Count')
plt.legend(loc=2)
plt.savefig('doc/figures/ufl_residual_final_distribution_10k.pdf')

# residual reduction
resreduce0 = sorted(np.loadtxt('data/resreduce0.txt', dtype=float))
resreduce1 = sorted(np.loadtxt('data/resreduce1.txt', dtype=float))
bin_edges   = np.logspace(-16, 0, 17)
bin_centers = np.linspace(-15.5, -0.5, 16)
hist0 = np.histogram(resreduce0, bin_edges)
hist1 = np.histogram(resreduce1, bin_edges)
cum0 = []
cum1 = []
rsum = 0
for a in hist0[0]: 
    rsum += a
    cum0.append(rsum)
rsum = 0
for a in hist1[0]: 
    rsum += a
    cum1.append(rsum)
print 'Residual reduction cumulation sum for USYMLQ = ', cum0
print 'Residual reduction cumulation sum for USYMQR = ', cum1
plt.figure()
ax = plt.subplot(111)
ax.bar(bin_centers+0.175, hist1[0], width=0.35, color='b', align='center', label='USYMQR')
ax.bar(bin_centers-0.175, hist0[0], width=0.35, color='r', align='center', label='USYMLQ')
ax.set_xticks(np.log10(bin_edges[1:]))
ax.set_xlabel('log(norm(r)/norm(b))')
ax.set_ylabel('Count')
plt.legend(loc=2)
plt.savefig('doc/figures/ufl_residual_reduction_distribution_10k.pdf')

# plt.show()
