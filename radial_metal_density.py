
import numpy as np
import matplotlib.pyplot as plt


names = ['dwSN', 'dwALL1', 'dwALL8']
galIDs = ['D9o2', 'D9q', 'D9m4a']
expns1 = ['0.900', '0.926', '0.950', '0.975', '0.990', '1.002']
expns2 = ['0.901', '0.925', '0.950', '0.976', '0.991', '1.001']
expns3 = ['0.900', '0.925', '0.950', '0.975', '0.990', '1.000']
expns = [expns1, expns2, expns3]
ions = ['HI', 'MgII', 'CIV', 'OVI']
atoms = ['H', 'Mg', 'C', 'O']

dstep = 0.1
dmin = [i*dstep for i in range(0,15)]
dmax = [i+dstep for i in dmin]
fileloc = '/home/jacob/research/dwarfs/gasfiles/'

for ion, atom  in zip(ions, atoms):

    fig1 = plt.figure(1)
    fig2 = plt.figure(2)

    print ion
    namecount = -1
    for galID, expn in zip(galIDs, expns):
        namecount+=1
        print galID
        r = []
        nIon = []
        nAtom = []
        for a in expn:

            filename = galID+'_GZa'+a+'.'+ion+'.txt'
            print '\t{0:s}'.format(filename)
        
            x, y, z, d, dense = np.loadtxt(fileloc+filename, 
                                               skiprows=2, 
                                               usecols=(1,2,3,11,13), 
                                               unpack=True)
            for i in range(0,len(x)):
                r.append(np.sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]))
                nIon.append(dense[i])
                nAtom.append(d[i])
            

        # Plot the mean densities
        nSum = np.zeros(len(dmin))
        nCount = np.zeros(len(dmin))
        nMean = np.zeros(len(dmin))
        rMean = np.zeros(len(dmin))
        nAtomSum = np.zeros(len(dmin))
        nAtomCount = np.zeros(len(dmin))
        nAtomMean = np.zeros(len(dmin))
        for i in range(0,len(r)):
            for j in range(0,len(dmin)):
                if r[i]>dmin[j] and r[i]<dmax[j]:
                    nSum[j] += nIon[i]
                    nCount[j] += 1.0
                    nAtomSum[j] += nAtom[j]
                    nAtomCount[j] += 1.0

        for i in range(0,len(nSum)):
            nMean[i] = np.log10(nSum[i] / nCount[i])
            nAtomMean[i] = np.log10(nAtomSum[i] / nAtomCount[i])

        for i in range(0,len(dmin)):
            rMean[i] = (dmin[i]+dmax[i])/2.0
        
        plt.figure(fig1.number) 
        plt.plot(rMean, nMean, marker='x', label=names[namecount])
        plt.figure(fig2.number)
        plt.plot(rMean, nAtomMean, marker='x', label=names[namecount])
    
    plt.figure(fig1.number)
    plt.xlabel('D/Rvir')
    plt.ylabel('Mean {0:s} Density [cm$^-3$]'.format(ion))
    plt.legend(frameon=False)
    plt.savefig('{0:s}_metalNumberDensity.png'.format(ion), bbox_inches='tight')
    plt.cla()
    plt.clf()

    plt.figure(fig2.number)
    plt.xlabel('D/Rvir')
    plt.ylabel('Mean {0:s} Density [cm$^-3$]'.format(atom))
    plt.legend(frameon=False)
    plt.savefig('{0:s}_metalNumberDensity.png'.format(atom), bbox_inches='tight')
    plt.cla()
    plt.clf()



