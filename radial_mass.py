#!/usr/bin/python

# Code to plot the radial gas mass profile for all three
# feedback mechanisms

import matplotlib.pyplot as plt
import numpy as np
import os
import sys



def convertDensity( numDense ):

    # Converts number density of an ion to mass density in 
    # solar masses per cubic parsec
    
    # Constants
    amu = 1.6605e-24  # in grams
    solar = 1.989e33  # in grams
    pc = 3.0856e18     # in cm

    atomicMass = 1.008
    
    # Unlog the number density
#    numDense = pow(10, numDense)

    # Convert to solar mass
    density = numDense * atomicMass * amu / solar

    # Convert to cubic parsec
    density = density * pow(pc,3)

    return density




galID_list = ['D9o2', 'D9q', 'D9m4a']
galdec_list = ['dwSN', 'dwALL_1', 'dwALL_8']
expn_list_1 = ['0.900', '0.926', '0.950', '0.975', '0.990', '1.002']
expn_list_2 = ['0.901', '0.925', '0.950', '0.976', '0.991', '1.001']
expn_list_3 = ['0.900', '0.925', '0.950', '0.975', '0.990', '1.000']
expn_list = [expn_list_1, expn_list_2, expn_list_3]
sym1 = ['x', 's', 'o']
sym2 = ['v', '^', 'D']
col = ['b', 'r', 'g']
rvir = [80.25497, 82.08840, 79.05807]
stellarMass = [8.1e8, 1.3e8, 2.1e7]
sfr = [0.07, 0.025, 0.0018]

minr = 0.0
maxr = 2.0
step = 0.1
Rmin = []
for i in np.arange(minr, maxr, step):
    Rmin.append(i)
Rmax = [i+step for i in Rmin]


fig,(p11, p21, p31) = plt.subplots(3,1,figsize=(5,10.8))
p12 = p11.twinx()
p22 = p21.twinx()
p32 = p31.twinx()
plotList = [p11, p21, p31]

fit, ax1 = plt.subplots()
ax2 = ax1.twinx()
for i in range(0,len(galID_list)):
    
    galID = galID_list[i]
    galdec = galdec_list[i]
    expn = expn_list[i]

    r = []
    n = []
    m = []
    size = []
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
        nH = data[:,7] 
        cellSize = data[:,0]
        for j in range(0,len(x)):
            dist = np.sqrt(x[j]*x[j] + y[j]*y[j] + z[j]*z[j])
            dist = dist / rvir[i]
            r.append(dist)
            
            # Calculate the mass in this cell
            mass = convertDensity(nH[j])*pow(cellSize[j], 3)
            m.append(mass)
            

    print 'Data Read In, Start Binning'
    # Bin the data
    print len(r)
    massmean, rmean, nerr = [],[],[]
    density = []
    for j in range(0,len(Rmin)):
        rmin = Rmin[j]
        rmax = Rmax[j]
        rsum, massSum = [], []
        for k in range(0,len(r)):
            if r[k]>rmin and r[k]<rmax:
                rsum.append(r[k])
                massSum.append(m[k])
        rmean.append(np.mean(rsum))
        binMass = np.mean(massSum)
        massmean.append(binMass)

        # Calcualte the volume of the sphere in this bin
        volume = 4/3 * np.pi * (pow(rmax*rvir[i],3)-pow(rmin*rvir[i],3))
        density.append(binMass/volume)

    print 'Binning Done, Begin Scaling'
    starScaleM, starScaleD = [], []
    sfrScaleM, sfrScaleD = [], []
    for j in range(0,len(rmean)):
        starScaleM.append( massmean[j] / stellarMass[i] )
        starScaleD.append( density[j] / stellarMass[i] )
        sfrScaleM.append( massmean[j] / sfr[i] )
        sfrScaleD.append( density[j] / sfr[i] )

    print 'Scaling Done, Begin Plotting'
    # Plot
    massmean = np.log10(massmean)
    density = np.log10(density)
    starScaleM = np.log10(starScaleM)
    starScaleD = np.log10(starScaleD)
    sfrScaleM = np.log10(sfrScaleM)
    sfrScaleD = np.log10(sfrScaleD)
    nerr = np.log10(nerr)
    nerr = [0.0 for s in nerr]
    print massmean
    print starScaleM
    print sfrScaleM 
    p11.plot( rmean, massmean, marker=sym1[i], color=col[i], ls='-', label=galdec+' mass')
    p12.plot( rmean, density, marker=sym2[i], color=col[i], ls='--', label=galdec+' density')
    
    p21.plot( rmean, starScaleM, marker=sym1[i], color=col[i], ls='-', label=galdec+' mass')
    p22.plot( rmean, starScaleD, marker=sym2[i], color=col[i], ls='--', label=galdec+' density')
    
    p31.plot( rmean, sfrScaleM, marker=sym1[i], color=col[i], ls='-', label=galdec+' mass')
    p32.plot( rmean, sfrScaleD, marker=sym2[i], color=col[i], ls='--', label=galdec+' density')
    


p11.legend(loc='upper right', frameon=False, fontsize='small')
p12.legend(loc='upper center', frameon=False, fontsize='small')
p21.legend(loc='upper right', frameon=False, fontsize='small')
p22.legend(loc='upper center', frameon=False, fontsize='small')
p31.legend(loc='upper right', frameon=False, fontsize='small')
p32.legend(loc='upper center', frameon=False, fontsize='small')

p11.set_xlabel('Distance [Rvir]')
p21.set_xlabel('Distance [Rvir]')
p31.set_xlabel('Distance [Rvir]')

#ax1.set_ylim([2, 4])
#ax1.set_ylim([-7,-3])

p21.set_ylabel('$\log($ $M_{gas}$ / $M_{*}$ )')
p22.set_ylabel('$\log($ $\rho_{gas}$ / $M_{*}$ )', labelpad=30, rotation=270)
 
p11.set_ylabel('$\log($ $M_{gas}$ [$M_{\odot}$]')
p12.set_ylabel('$\log($ $\rho_{gas}$ [$M_{\odot}$/kpc$^{3}$] )', labelpad=30, rotation=270)

p31.set_ylabel('$\log($ $M_{gas}$ / <SFR> )')
p32.set_ylabel('$\log($ $\rho_{gas}$ / <SFR> )', labelpad=30, rotation=270)
#fig.tight_layout()
fig.savefig('master_radial_mass_bulk.pdf', bbox_inches='tight')
    

