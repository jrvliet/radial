# 
# Plots the stellar raidal profile
#
import sys
import yt
from yt.fields.api import ValidateParameter
import numpy as np
import matplotlib.pyplot as plt



def _gc_distance(field, data):
    # A new field for the distance
    # between the gas cell and the 
    # location of maximum density
    center = data.get_field_parameter('center')
    xdist = data[('stars', 'particle_position_x')] - center[0]
    ydist = data[('stars', 'particle_position_y')] - center[1]
    zdist = data[('stars', 'particle_position_z')] - center[2]
    r = np.sqrt( xdist*xdist + ydist*ydist + zdist*zdist )
    return r




fileLoc = '/praesepe/jrvander/dataFiles/'
fileLoc = '/fomalhaut-data/jrvander/simulations/dwarfs/dataFiles/'

galID_list = ['D9o2', 'D9q', 'D9m4a']
labels = ['dwSN', 'dwALL\_1', 'dwALL\_8']
col = ['k', 'b', 'r']
sym = ['s', '^', 'o']
expn_list_1 = ['0.900', '0.926', '0.950', '0.975', '0.990', '1.002']
expn_list_2 = ['0.901', '0.925', '0.950', '0.976', '0.991', '1.001']
expn_list_3 = ['0.900', '0.925', '0.950', '0.975', '0.990', '1.000']
expn_list = [expn_list_1, expn_list_2, expn_list_3]
rvir = [80.25497, 82.08840, 79.05807]
rvir = np.asarray(rvir)

minr = 0.0
maxr = 2.0
step = 0.1
bins = np.arange(minr, maxr, step)

fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()

outputFile = 'binnedStars.out'
outf = open(fileLoc+outputFile, 'w')

for i in range(0,len(galID_list)):
    
    galID = galID_list[i]
    expn = expn_list[i]

    allDist = np.zeros(20)
    print galID    
    
    for a in expn:
        print a
        galFile = fileLoc+galID+'/10MpcBox_HartGal_csf_a'+a+'.d'

        ds = yt.load(galFile)

        ds.add_field( 'gc_distance', 
              function = _gc_distance, 
              units = 'cm',
              take_log = False,
              validators = [ValidateParameter('center') ] )

        # Index the data structure
        ad = ds.all_data()

        # Create a sphere centered around the galaxy of radius 2*Rvir
        sp = ds.sphere('max', (2*rvir[i], 'kpc'))
    
        r = sp[('gc_distance')].in_units('kpc')
        r = r/rvir[i]
    
        # Bin the data
        inds = np.digitize(r, bins)
        dist = np.bincount(inds)
        if len(dist)>len(allDist): 
            allDist = np.add( allDist, dist[:-1] )
        else: 
            allDist = np.add(allDist, dist)

        ad.clear_data()

    allDist = allDist / 6.    # Take the average
    allDist = np.log10( allDist )
    plt.plot(bins, allDist, marker=sym[i], color=col[i], label=labels[i])

plt.xlabel('Distance [Rvir]')
plt.ylabel('Log ( Number of Stars )')
plt.legend(frameon=False)
plt.savefig('stellarProfile.pdf')



