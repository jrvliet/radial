#!/usr/bin/python

# Code to plot the radial temperature profile for all three
# feedback mechanisms

import matplotlib.pyplot as plt
import numpy as np
import os
import sys

#record = int(sys.argv[1])
record = 0

galID_list = ['D9o2', 'D9q', 'D9m4a']
galdec_list = ['dwSN', 'dwALL_1', 'dwALL_8']
expn_list_1 = ['0.900', '0.926', '0.950', '0.975', '0.990', '1.002']
expn_list_2 = ['0.901', '0.925', '0.950', '0.976', '0.991', '1.001']
expn_list_3 = ['0.900', '0.925', '0.950', '0.975', '0.990', '1.000']
expn_list = [expn_list_1, expn_list_2, expn_list_3]
sym = ['x-', 's-', 'o-']
rvir = [80.25497, 82.08840, 79.05807]
mass = [8.1e8, 1.3e8, 2.1e7]

minr = 0.0
maxr = 2.0
step = 0.1
Rmin = []
for i in np.arange(minr, maxr, step):
    Rmin.append(i)
Rmax = [i+step for i in Rmin]

lines = []
radius = []



meanfig = plt.figure()
meanax = meanfig.add_subplot(111)
medfig = plt.figure()
medax = medfig.add_subplot(111)


for i in range(0,len(galID_list)):
    
    galID = galID_list[i]
    galdec = galdec_list[i]
    expn = expn_list[i]

    r = []
#    n = []
    T = []
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
#        density = data[:,7]
        temp = data[:,8]

        for j in range(0,len(x)):
            
            dist = np.sqrt(x[j]*x[j] + y[j]*y[j] + z[j]*z[j])
            dist = dist / rvir[i]
            r.append(dist)
#            n.append(density[j])
            T.append(temp[j])

        # Bin the data
        rmean, Tmean, Terr = [],[],[]
        rmed, Tmed = [], []
        for j in range(0,len(Rmin)):
        
            rmin = Rmin[j]
            rmax = Rmax[j]
            rsum, Tsum = [], []
            for k in range(0,len(r)):
                if r[k]>rmin and r[k]<rmax:
                    rsum.append(r[k])
                    Tsum.append(T[k])

            rmed.append(np.median(rsum))
            Tmed.append(np.median(Tsum))
            rmean.append(np.mean(rsum))
            Tmean.append(np.mean(Tsum))
            Terr.append(np.std(Tsum))

        lines.append( Tmean )
        radius = np.copy( rmean )

#        if i==0:
#            scale = np.copy(Tmean)
#        for j in range(0,len(Tmean)):
#            Tmean[j] = Tmean[j]/scale[j]


                   
        # Plot
        Tmean = np.log10(Tmean)
        Tmed = np.log10(Tmed)
        Terr = np.log10(Terr)
        Terr = [0.0 for s in Terr]
        meanax.errorbar(rmean, Tmean, yerr=Terr, fmt=sym[i], label=galdec)
        medax.errorbar(rmed, Tmed, yerr=Terr, fmt=sym[i], label=galdec)
        

#plt.legend(loc='lower right', frameon=False)
medax.set_xlabel('Distance [Rvir]')
medax.set_ylabel('log (T [K])')
medfig.savefig('master_radial_temp_median.pdf', bbox_inches='tight')

meanax.set_xlabel('Distance [Rvir]')
meanax.set_ylabel('log (T [K])')
meanfig.savefig('master_radial_temp_median.pdf', bbox_inches='tight')





if record==1:
    print 'Writing to file'
    # Write to file
    outputF = 'binnedTemp.out'
    outf = open(outputF, 'w')
    outf.write('R \t D9o2 \t D9q \t D9m4a \n')
    for i in range(0,len(lines[0])):

        s = '{0:.1f} \t {1:.4f} \t {2:.4f} \t {3:.4f} \n'.format(radius[i], np.log10(lines[0][i]), np.log10(lines[1][i]), np.log10(lines[2][i]))
        outf.write(s)

    outf.close()





