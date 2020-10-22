#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import glob

mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
#mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

sta = 'BFO'
fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize=(12,12))
ax = [ax1, ax2, ax3]
for idx2, tidetype in enumerate(['_semidiurnal', '_diurnal']):
    filename = glob.glob('*' + sta + '*' + tidetype)[0]
    print(filename)
    results = []
    with open(filename) as f:
        next(f)
        for line in f:
            line = line.rstrip()
            line = line.replace(' ','')
            line = line.split(',')
            results.append(line)


    vals = {}

    for idx, comp in enumerate(['Z', 'N', 'E']):
        rats, cols, times = [], [], []
        for result in results:
            if result[0] == comp:
                rat = float(result[3])/float(result[4])
                if rat > 1.3:
                    rat = 1.3
                if rat < 0.7:
                    rat = 0.7
                rats.append(rat)
                cols.append(float(result[2]))
                times.append(int(result[5]))
        if idx2 == 0:
            m = 'o'
        else:
            m ='s'
        ax[idx].scatter(times,rats, marker=m, c = cols, label= tidetype.replace('_','') + ' ' + comp + 
            ' Ratio: ' + str(round(np.mean(rats),3)), vmin=0.9, vmax=1.1)
        ax[idx].set_ylim((0.7,1.3))
        ax[idx].legend(loc=2)
        if idx == 1:
            ax[idx].set_ylabel('Amplitude Ratio')
        #sys.exit()
ax[2].set_xlabel('Time (DOY)')
ax[0].set_title(sta)
fig.savefig('Test.png', format='PNG')