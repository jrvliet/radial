#!/usr/bin/python

# Code to plot the radial temperature profile for all three
# feedback mechanisms

import matplotlib.pyplot as plt
import numpy as np
import os
import sys

scaled = 0     # Determines if mass fraction is scaled by stellar mass

galID_list = ['D9o2', 'D9q', 'D9m4a']
galdec_list = ['dwSN', 'dwALL_1', 'dwALL_8']
expn_list_1 = ['0.900', '0.926', '0.950', '0.975', '0.990', '1.002']
expn_list_2 = ['0.901', '0.925', '0.950', '0.976', '0.991', '1.001']
expn_list_3 = ['0.900', '0.925', '0.950', '0.975', '0.990', '1.000']
expn_list = [expn_list_1, expn_list_2, expn_list_3]
sym1 = ['s', '^', 'o']
sym2 = ['s', '^', 'o']
col = ['k', 'b', 'r']    # Plotting colors
rvir = [80.25497, 82.08840, 79.05807]
mass = [8.1e8, 1.3e8, 2.1e7]
sfr = [0.07, 0.025, 0.0018]
minr = 0.0
maxr = 2.0
step = 0.1
Rmin = []
for i in np.arange(minr, maxr, step):
    Rmin.append(i)
Rmax = [i+step for i in Rmin]

fig,(p1, p2, p3) = plt.subplots(3,1,figsize=(6.5,10.8))
plotList = [p1, p2, p3]
p12 = p1.twinx()
p22 = p2.twinx()
p32 = p3.twinx()
#fit, ax1 = plt.subplots()
#ax2 = ax1.twinx()
for i in range(0,len(galID_list)):
    
    galID = galID_list[i]
    galdec = galdec_list[i]
    expn = expn_list[i]

    r = []
    n = []
    snII = []
    snIa = []
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
        sn2 = data[:,9] 
        sn1 = data[:,10]

        for j in range(0,len(x)):
            
            dist = np.sqrt(x[j]*x[j] + y[j]*y[j] + z[j]*z[j])
            dist = dist / rvir[i]
            r.append(dist)
            snII.append(sn2[j])
            snIa.append(sn1[j])

    print 'Data Read In, Start Binning'
    # Bin the data
    print len(r)
    snIamean, snIImean, rmean, nerr = [],[],[],[]
    for j in range(0,len(Rmin)):
        
        rmin = Rmin[j]
        rmax = Rmax[j]
        rsum, snIaSum, snIISum = [], [], []
        for k in range(0,len(r)):
            if r[k]>rmin and r[k]<rmax:
                rsum.append(r[k])
                snIaSum.append(snIa[k])
                snIISum.append(snII[k])
        rmean.append(np.mean(rsum))
        snIImean.append(np.mean(snIISum))
        snIamean.append(np.mean(snIaSum))

    
    print 'Binning Done, Begin Scaling'
    starScaleII, starScaleIa = [], []
    sfrScaleII, sfrScaleIa = [], []
    for j in range(0,len(rmean)):
        starScaleII.append( snIImean[j] / mass[i] )
        starScaleIa.append( snIamean[j] / mass[i] )
        sfrScaleII.append( snIImean[j] / sfr[i] )
        sfrScaleIa.append( snIamean[j] / sfr[i] )
    print 'Scaling Done, Begin Plotting'
    # Plot
    snIImean = np.log10(snIImean)
    snIamean = np.log10(snIamean)

    starScaleII = np.log10(starScaleII)
    starScaleIa = np.log10(starScaleIa)

    sfrScaleII = np.log10(sfrScaleII)
    sfrScaleIa = np.log10(sfrScaleIa)

    nerr = np.log10(nerr)
    nerr = [0.0 for s in nerr]

#    ax1 = plotList[i]
#    ax2 = ax1.twinx()
#    ax2.plot( rmean, snIamean, marker=sym2[i], color=col[i], ls='--', label=galdec+', SnIa')
    p1.plot( rmean, snIImean, marker=sym1[i], color=col[i], ls='-', label=galdec+', SnII')
    p12.plot( rmean, snIamean, marker=sym2[i], color=col[i], ls='--', label=galdec+', SnIa', mfc='none', mec=col[i])

    p2.plot( rmean, starScaleII, marker=sym1[i], color=col[i], ls='-', label=galdec+', SnII')
    p22.plot( rmean, starScaleIa, marker=sym2[i], color=col[i], ls='--', label=galdec+', SnIa', mfc='none', mec=col[i])

    p3.plot( rmean, sfrScaleII, marker=sym1[i], color=col[i], ls='-', label=galdec+', SnII')
    p32.plot( rmean, sfrScaleIa, marker=sym2[i], color=col[i], ls='--', label=galdec+', SnIa', mfc='none', mec=col[i])



p1.legend(loc='upper right', frameon=False, fontsize='small')
p12.legend(loc='upper center', frameon=False, fontsize='small')
p2.legend(loc='upper right', frameon=False, fontsize='small')
p22.legend(loc='upper center', frameon=False, fontsize='small')
p3.legend(loc='upper right', frameon=False, fontsize='small')
p32.legend(loc='upper center', frameon=False, fontsize='small')


p1.set_xlabel('Distance [Rvir]')
p2.set_xlabel('Distance [Rvir]')
p3.set_xlabel('Distance [Rvir]')

p3.set_ylim([-2.5, 0.0])
p32.set_ylim([-4.0, -1.5])

p1.set_ylabel('SNII MF')
p12.set_ylabel('SNIa MF', labelpad=30, rotation=270)

p2.set_ylabel('SNII MF / $M_{*}$')
p22.set_ylabel('SNIa MF / $M_{*}$', labelpad=30, rotation=270)

p3.set_ylabel('SNII MF / <SFR>')
p32.set_ylabel('SNIa MF / <SFR>', labelpad=30, rotation=270)

fig.tight_layout()
fig.savefig('master_radial_metal_bulk.pdf', bbox_inches='tight')
    

