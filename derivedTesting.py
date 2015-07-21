import yt
from yt.fields.api import ValidateParameter
import numpy as np
import matplotlib.pyplot as plt

def _gc_distance(field, data):

#    if data.has_field_parameter('center'):
#        center = data.get_field_parameter("center").in_units("cm")
#    else:
#        center = data.ds.arr(np.zeros(3), "cm")

    center = data.get_field_parameter('center')
    xdist = data[('stars', 'particle_position_x')] - center[0]
    ydist = data[('stars', 'particle_position_y')] - center[1]
    zdist = data[('stars', 'particle_position_z')] - center[2]
    r = np.sqrt( xdist*xdist + ydist*ydist + zdist*zdist )
    return r


fileLoc = '/fomalhaut-data/jrvander/simulations/dwarfs/dataFiles/'
rvir = 79.05807
galFile = fileLoc+'D9m4a/10MpcBox_HartGal_csf_a1.000.d'
ds = yt.load(galFile)
ds.add_field( 'gc_distance', function=_gc_distance, units='cm', take_log=False, validators=[ValidateParameter('center')] )
ds.index
ad = ds.all_data()
sp = ds.sphere('max', (2*rvir, "kpc"))

r = sp[('gc_distance')].in_units('kpc')

r = r/rvir


minr = 0.0
maxr = 2.0
step = 0.1

bins = np.arange(minr, maxr, step)
inds = np.digitize(r, bins)
dist = np.bincount(inds)
dist = np.log10( dist )
#plt.hist(r, bins=20, histtype='step')
print bins
print len(bins)
print dist
print len(dist)
plt.plot(bins, dist, 'rx')
plt.savefig('compare.pdf')

ad.clear_data()

