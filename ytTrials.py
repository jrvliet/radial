# coding: utf-8
import yt
testLoc = '/home/jacob/matrix/sebass_gals/dwarfs/testGalaxy/'
galFile = testLoc+'D9p_500/10MpcBox_HartGal_csf_a0.500.d'
ds = yt.load(galFile)
ad = ds.all_data()
center = ds.find_max(('gas', 'density'))
print center
centerVal, centerLoc = ds.find_max(('gas', 'density'))
print centerLoc
centerVal, centerLoc = ds.find_max(('gas', 'density')).in_units("kpc")
centerLoc.in_units('kpc')
sp = ds.sphere(centerLoc, (85, "kpc"))
x = ad[("stars", "particle_position_x")].in_units("kpc")
x = sp[("stars", "particle_position_x")].in_units("kpc")
xSP = sp[("stars", "particle_position_x")].in_units("kpc")
x = ad[("stars", "particle_position_x")].in_units("kpc")
print len(x)
print len(xSP)
print x
if xSP not in x:
    print xSP
    
for i in xSP:
    if i not in x:
        print i
        
for i in xSP:
    if i not in x:
        print i
        
for i in x:
    if i not in xSP:
        print i
        
for i in x:
    if i not in xSP:
        print i-6754.03795161
        
for i in x:
    if i not in xSP:
        print i-centerLoc[0]
        
x = ad[("stars", "particle_position_x")].in_units("kpc")
y = ad[("stars", "particle_position_y")].in_units("kpc")
z = ad[("stars", "particle_position_z")].in_units("kpc")
r = []
from numpy import sqrt
for i in x:
    for j in y:
        for k in z:
            dist = sqrt(i*i + j*j + k*k)
            r.append(dist)
            
r = []
for i in range(0,len(x)):
    r.append( sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i] ))
    
print max(r)
print min(r)
x = [i-centerLoc[0] for i in x]
y = [i-centerLoc[1] for i in y]
z = [i-centerLoc[2] for i in z]
for i in range(0,len(x)):
    r.append( sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i] ))
    
print min(r)
print max(r)
r = []
for i in range(0,len(x)):
    r.append( sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i] ))
    
print max(r)
sx = [i-centerLoc[0] for i in xSP]
sy = [i-centerLoc[1] for i in ySP]
xSP = sp[("stars", "particle_position_x")].in_units("kpc")
sy = sp[("stars", "particle_position_y")].in_units("kpc")
sz = sp[("stars", "particle_position_z")].in_units("kpc")
sy = [i-centerLoc[1] for i in sy]
sz = [i-centerLoc[2] for i in sz]
sr = []
for i in range(0,len(sx)):
    sr.append( sqrt(sx[i]*sx[i] + sy[i]*sy[i] + sz[i]*sz[i] ))
    
print min(sr)
print max(sr)
import matplotlib.pyplot as plt
plt.plot(r)
plot.plot(sr)
plt.plot(sr)
plt.savefig('dum.png')
yt.ProjectionPlot(ds, "x", "temperature", weight_field="density").save()
yt.ProjectionPlot(ds, "x", "density", weight_field="density").save()
yt.ProjectionPlot(sp, "x", "density", weight_field="density").save()
plot = yt.ProfilePlot(my_sphere, "radius", "density")
plot = yt.ProfilePlot(sp, sr, ('stars','particle_mass'))
plt.plot(sr, '.')
plt.cla()
plt.clf()
plt.plot(sr, '.')
plt.savefig('dum2.png')
plt.hist(sr, bins=20)
plt.cla()
plt.clf()
plt.hist(sr, bins=20)
plt.savefig('dumhist.png')
plt.savefig('dumhist.png')
plot = yt.ProjectionPlot(sp, 'x', ('stars', 'particle_mass'))
leftCorner = centerLoc - ds.arr([100, 100, 100], 'kpc')
rightCorner = centerLoc + ds.arr([100, 100, 100], 'kpc')
region = ds.box(leftCorner, rightCorner)
plot = yt.ProjectionPlot(region, 'x', ('stars', 'particle_mass'))
ds.datasource = region
p = yt.ProjectionPlot(ds,'x',[('stars','particle_mass')], data_source=region, center=center, width=ds.arr([1,1],'Mpc'))
p = yt.ProjectionPlot(ds,'x',('stars','particle_mass'), data_source=region, center=center, width=ds.arr([1,1],'Mpc'))
p = yt.ProjectionPlot(ds,'x',('stars','particle_mass'), data_source=region, center=centerLoc, width=ds.arr([1,1],'Mpc'))
p = yt.ProjectionPlot(ds,'x',[('stars','particle_mass')], data_source=region, center=centerLoc, width=ds.arr([1,1],'Mpc'))
p = yt.ProjectionPlot(ds,'x',[('stars','particle_mass')], data_source=region, center=centerLoc, width=ds.arr([100,100],'kpc'))
p = yt.ProjectionPlot(ds,'x',[('stars','particle_mass')], data_source=region, center=centerLoc, width=(100, 'kpc'))
p = yt.ProjectionPlot(ds,'x',[('stars','particle_mass')], data_source=region, center=center, width=(100, 'kpc'))
p = yt.ProjectionPlot(ds, 'x', [('stars','particle_mass')], data_source=region, center=center, width=(60, 'kpc'))
p = yt.ProjectionPlot(ds, 'x', [('stars','particle_mass')], data_source=region, center='max', width=(60, 'kpc'))
p = yt.ProjectionPlot(ds, 'x', [('stars','particle_mass')], data_source=region, center='max', width=(100, 'kpc'))
p = yt.ProjectionPlot(ds, 'x', [('gas','density')], data_source=region, center='max', width=(100, 'kpc'))
p.save()
p.annotate_particles()
p.save()
p.save(1)
p.save('dum.png')
p.annotate_particles(1)
p.save('dum.png')
p1 = yt.ProjectionPlot(ds, 'x', [('gas','density')], data_source=region, center='max', width=(100, 'kpc'))
p1.annotate_particles()
p1.save('dum.png')
p1.annotate_particles((100, 'kpc'))
p1.save('dum.png')
p1.save()
p1.save('annotate.png', bbox_inches=’tight’)
p1.save('annotate.png', 'bbox_inches':’tight’)
p1.save('annotate.png', 'bbox_inches'=’tight’)
p1.save('annotate.png', {'bbox_inches':’tight’})
p1.save('annotate.png', mpl_kwargs={'bbox_inches':’tight’})
p1.save('annotate.png', mpl_kwargs={'bbox_inches':'tight’})
p1.save('annotate.png', 'bbox_inches':'tight’})
p1.save('annotate.png', 'bbox_inches'='tight’})
p1.save('annotate.png', bbox_inches='tight’})
p1.save('annotate.png', bbox_inches='tight’)
p1.save('annotate.png', bbox_inches='tight')
p1.save('annotate.png')
p1.save(slef, 'annotate.png')
p1.save(self, 'annotate.png')
p2 = yt.ProjectionPlot(ds, 'x', [('gas','density')], center='max', width=(100, 'kpc'))
p2.save()
p2.annotate()
p2.annotate_particles()
p2.save()
p2.save()
p3 = yt.ProjectionPlot(ds, 'x', [('gas','density')])
p3.annotate_particles()
p3.save()
ds.derived_field_list
#Define temperature derived field
def _temperature(field, data):
            r0 = data.ds.parameters['boxh'] / data.ds.parameters['ng']
            T0 = 3.03e5 * r0**2 * data.ds.parameters['wmu'] * data.ds.parameters['Om0']
            T0 = T0 * (data.ds.parameters['gamma']-1.) / (data.ds.parameters['aexpn']**2)
            T_conv = data.ds.quan(T0, 'K/code_velocity**2')
            return T_conv * data['art', 'GasEnergy'].d / data['art', 'Density'].d
ds.add_field(('gas', 'temperature'), function=_temperature, units='K')
ds.derived_field_list
p4 = yt.ProjectionPlot(ds, 'x', [('gas','temperature')])
def _temperature(field, data):
            r0 = data.ds.parameters['boxh'] / data.ds.parameters['ng']
            T0 = 3.03e5 * r0**2 * data.ds.parameters['wmu'] * data.ds.parameters['Om0']
            T0 = T0 * (data.ds.parameters['gamma']-1.) / (data.ds.parameters['aexpn']**2)
            T_conv = data.ds.quan(T0, 'K/code_velocity**2')
            return T_conv * data['art', 'GasEnergy'] / data['art', 'Density']
ds.add_field(('gas', 'temperature'), function=_temperature, units='K')
ds.add_field(('gas', 'temperature'), function=_temperature, units='K', force_override=True)
p4 = yt.ProjectionPlot(ds, 'x', [('gas','temperature')])
p4.save('temp.png')
get_ipython().magic(u'save')
get_ipython().magic(u'save ytTrials')
get_ipython().magic(u'save ytTrials.txt')
get_ipython().magic(u"save 'ytTrials.txt'")
get_ipython().magic(u'save ytTrials 1-145')
