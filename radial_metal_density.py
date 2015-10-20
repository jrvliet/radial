
import numpy as np
import matplotlib.pyplot as plt


names = ['dwSN', 'dwALL1', 'dwALL8']
galIDs = ['D9o2', 'D9q', 'D9m4a']
expns1 = ['0.900', '0.926', '0.950', '0.975', '0.990', '1.002']
expns2 = ['0.901', '0.925', '0.950', '0.976', '0.991', '1.001']
expns3 = ['0.900', '0.925', '0.950', '0.975', '0.990', '1.000']
expns = [expns1, expns2, expns3]
ions = ['HI', 'MgII', 'CIV', 'OVI']
ions = ['OVI']

dstep = 0.1
dmin = [i*dstep for i in range(0,15)]
dmax = [i+dstep for i in dmin]
fileloc = '/home/jacob/research/dwarfs/gasfiles/'

for ion in ions:
    print ion
    namecount = -1
    for galID, expn in zip(galIDs, expns):
        namecount+=1
        print galID
        r = []
        nIon = []
        for a in expn:

            filename = galID+'_GZa'+a+'.'+ion+'.txt'
            print '\t{0:s}'.format(filename)
        
            x, y, z, dense = np.loadtxt(fileloc+filename, skiprows=2, usecols=(1,2,3,13), unpack=True)
            print dense[0] 
            for i in range(0,len(x)):
                r.append(np.sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]))
                nIon.append(dense[i])
            

        # Plot the mean densities
        nSum = np.zeros(len(dmin))
        nCount = np.zeros(len(dmin))
        nMean = np.zeros(len(dmin))
        rMean = np.zeros(len(dmin))
        for i in range(0,len(r)):
            for j in range(0,len(dmin)):
                if r[i]>dmin[j] and r[i]<dmax[j]:
                    nSum[j] += nIon[i]
                    nCount[j] += 1.0

        for i in range(0,len(nSum)):
            nMean[i] = np.log10(nSum[i] / nCount[i])

        for i in range(0,len(dmin)):
            rMean[i] = (dmin[i]+dmax[i])/2.0
    
        plt.plot(rMean, nMean, marker='x', label=names[namecount])

    plt.xlabel('D/Rvir')
    plt.ylabel('Mean {0:s} Density [cm$^-3$]'.format(ion))
    plt.legend(frameon=False)
    plt.savefig('{0:s}_metalNumberDensity.png'.format(ion), bbox_inches='tight')
    plt.cla()
    plt.clf()
    



