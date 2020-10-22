#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import glob

mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
#mpl.rc('text', usetex=True)
mpl.rc('font',size=16)

tidetype = '_semidiurnal'
fig, ax = plt.subplots(1,3, figsize=(14,10))
ax = ax.flatten()
stas = []
files = glob.glob('IU*' + tidetype)
files.sort()
ratsall = []
for sidx, filename in enumerate(files):
    results = []
    with open(filename) as f:
        next(f)
        for line in f:
            line = line.rstrip()
            line = line.replace(' ','')
            line = line.split(',')
            results.append(line)

    vals = {}
    stas.append(filename.split('_')[1])
    for idx, comp in enumerate(['Z', 'N', 'E']):
        rats, cols, times = [], [], []
        for result in results:
            if result[0] == comp:
                rat = float(result[3])/float(result[4])
                if rat > 1.5:
                    rat = 1.5
                if rat < 0.5:
                    rat = 0.5
                if float(result[2]) > 0.97:
                    rats.append(rat)
                    ratsall.append(rat)

        rats = np.array(rats)
        ax[idx].errorbar(np.median(rats),sidx+0.5,xerr=np.std(rats), color ='k', marker='o')
        #ax[idx].set_ylim((0.7,1.3))
        #ax[idx].legend(loc=2)
        #if idx == 1:
        #    ax[idx].set_ylabel('Amplitude Ratio')
        #sys.exit()
print(stas)
print(str(np.mean(ratsall)) + ' ' + str(np.std(ratsall)))
letters = ['(a)', '(b)', '(c)']
for idx in range(3):
    ax[idx].set_yticks(np.array(range(len(stas)))+0.5)
    ax[idx].set_yticklabels(stas)
    ax[idx].set_xlim((0.5,1.5))
    ax[idx].set_ylim((0.5, len(stas)+0.5))
    ax[idx].text(0.3, len(stas)+1, letters[idx] )
ax[1].set_xlabel('Amplitude Ratio')
ax[0].set_ylabel('Station Code')
plt.tight_layout()
#ax[2].set_xlabel('Time (DOY)')
#ax[0].set_title(sta)
fig.savefig('Test.png', format='PNG')