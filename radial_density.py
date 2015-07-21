#!/usr/bin/python

# Code to plot the radial temperature profile for all three
# feedback mechanisms

import matplotlib.pyplot as plt
import numpy as np
import os
import sys

record = int(sys.argv[1])


galID_list = ['D9o2', 'D9q', 'D9m4a']
galdec_list = ['dwSN', 'dwALL_1', 'dwALL_8']
expn_list_1 = ['0.900', '0.926', '0.950', '0.975', '0.990', '1.002']
expn_list_2 = ['0.901', '0.925', '0.950', '0.976', '0.991', '1.001']
expn_list_3 = ['0.900', '0.925', '0.950', '0.975', '0.990', '1.000']
expn_list = [expn_list_1, expn_list_2, expn_list_3]
sym = ['x-', 's-', 'o-']
rvir = [80.25497, 82.08840, 79.05807]
minr = 0.0
maxr = 2.0
step = 0.1
Rmin = []
for i in np.arange(minr, maxr, step):
    Rmin.append(i)
Rmax = [i+step for i in Rmin]

lines = []
radius = []

for i in range(0,len(galID_list)):
    
    galID = galID_list[i]
    galdec = galdec_list[i]
    expn = expn_list[i]

    r = []
    n = []
#    T = []
    print ''
    print galID
    for a in expn:
        print a
        filepath = '/home/matrix3/jrvander/sebass_gals/dwarfs/'+galID+'_outputs/a'+a+'/'

        filename = galID+'_GZa'+a+'.txt'

        data = np.loadtxt(filepath+filename, skiprows=2)
        x = data[:,1]
        y = data[:,2]
        z = data[:,3]
        density = data[:,7]

        for j in range(0,len(x)):
            
            dist = np.sqrt(x[j]*x[j] + y[j]*y[j] + z[j]*z[j])
            dist = dist / rvir[i]
            r.append(dist)
            n.append(density[j])

        # Bin the data
        rmean, nmean, nerr = [],[],[]
        for j in range(0,len(Rmin)):
        
            rmin = Rmin[j]
            rmax = Rmax[j]
            rsum, nsum = [], []
            for k in range(0,len(r)):
                if r[k]>rmin and r[k]<rmax:
                    rsum.append(r[k])
                    nsum.append(n[k])

            rmean.append(np.mean(rsum))
            nmean.append(np.mean(nsum))
            nerr.append(np.std(nsum))


        lines.append( nmean )
        radius = np.copy( rmean )
    
        # Plot
        if i==0:
            scale = np.copy(nmean)
        for j in range(0,len(nmean)):
            nmean[j] = nmean[j]/scale[j]
        nmean = np.log10(nmean)
        nerr = np.log10(nerr)
        nerr = [0.0 for s in nerr]
        plt.errorbar(rmean, nmean, yerr=nerr, fmt=sym[i], label=galdec)
    
#    plt.legend(loc='upper right', frameon=False)
    plt.xlabel('Distance [Rvir]')
    plt.ylabel('log (n$_{H}$ [cm$^{-3}$])')
    if record==1:
        plt.savefig('master_radial_density.eps', bbox_inches='tight')
    else:
        plt.savefig('master_radial_density_bulk.pdf', bbox_inches='tight')

if record==1:
    print 'Writing to file'
    # Write to file
    outputF = 'binnedDensity.out'
    outf = open(outputF, 'w')
    outf.write('R \t D9o2 \t D9q \t D9m4a \n')
    for i in range(0,len(lines[0])):

        s = '{0:.1f} \t {1:.4f} \t {2:.4f} \t {3:.4f} \n'.format(radius[i], np.log10(lines[0][i]), np.log10(lines[1][i]), np.log10(lines[2][i]))
        outf.write(s)

    outf.close()







