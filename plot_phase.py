#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import glob

mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
#mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

sta = 'CCM'
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
        phases, cols, times = [], [], []
        for result in results:
            if result[0] == comp:
                phase = float(result[6])*180./np.pi
                phases.append(phase)
                cols.append(float(result[2]))
                times.append(int(result[8]))
        if idx2 == 0:
            m = 'o'
        else:
            m ='s'
        ax[idx].scatter(times,phases, marker=m, c = cols, label= tidetype.replace('_','') + ' ' + comp + 
            ' Phase: ' + str(round(np.mean(phases),3)), vmin=0.9, vmax=1.1)
        #ax[idx].set_ylim((0.7,1.3))
        ax[idx].legend(loc=2)
        if idx == 1:
            ax[idx].set_ylabel('Phase Difference (degrees)')
        #sys.exit()
ax[2].set_xlabel('Time (DOY)')
ax[0].set_title(sta)
fig.savefig('Test_Phase.png', format='PNG')