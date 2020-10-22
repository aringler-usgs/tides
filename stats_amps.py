#!/usr/bin/env python
import glob
import sys
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import matplotlib.font_manager
import cartopy as cart

# Importing and applying font
mpl.rc('font', family = 'serif')
mpl.rc('font', serif = 'Times') 
mpl.rc('text', usetex = True)
mpl.rc('font', size=18)
net = 'IU'
fig = plt.figure(1, figsize=(14,12))
stime = UTCDateTime('2019-001T00:00:00')
etime = UTCDateTime('2019-010T00:00:00')
client = Client('IRIS')

for idx2, tidetype in enumerate(['semidiurnal', 'diurnal']):
    stas = glob.glob(net + '*_' + tidetype)
    stas = [sta.split('_')[1] for sta in stas]
    amps_good, amps_goodstd, lens_good = {}, {}, {}

    for sta in stas:
        filename = glob.glob(net + '*' + sta + '*' + tidetype)[0]
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

            amps, cols, times = [], [], []
            for result in results:
                if (result[0] == comp) and (float(result[2]) > .95):
                    amp = float(result[3])/float(result[4])
                    #phase = float(result[6])*180/np.pi
                    if amp > 1.5:
                        amp = 1.5
                    if amp < .5:
                        amp = 0.5
                    amps.append(amp)
            lens_good[sta + '_' + comp] = len(amps)
            amps_good[sta + '_' + comp] = np.mean(amps)
            amps_goodstd[sta + '_' + comp] = np.std(amps)


    letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
    for idx, comp in enumerate(['Z', 'N', 'E']):
        vals, std, lens = [], [], []
        for sta in stas:
            try:
                vals.append(amps_good[sta + '_' + comp])
                std.append(amps_goodstd[sta + '_' + comp])

                lens.append(100*lens_good[sta + '_' + comp]/365.)
            except:
                continue

        ax = fig.add_subplot(3,2,2*idx+idx2+ 1)
        ax.hist(vals, bins=10)
        ax.set_ylim((0.,40.))
        ax.text(.25, 40, letters[2*idx + idx2 ] )

ax = fig.add_subplot(3,2,3)
ax.set_ylabel('Hits')
ax = fig.add_subplot(3,2,5)
ax.set_xlabel('Amplitude Ratio')
ax = fig.add_subplot(3,2,6)
ax.set_xlabel('Amplitude Ratio')

        #ax.coastlines()
        #ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
        #ax.set_global()
        #ax.coastlines()
        #
        #im = ax.scatter(alons, alats, c=vals, transform=ccrs.Geodetic(), zorder=3, vmin=-10, vmax = 10)
        

#cb_ax = fig.add_axes([0.2, 0.05, 0.6, 0.03])
#cbar = fig.colorbar(im, cax=cb_ax, orientation='horizontal')
#cbar.set_label('Phase Difference (degree)') 
#pos = cbar.ax.get_position()
#ax2 = cbar.ax.twiny()
#ax2.set_xlim((-10.,10.))
#ax2.set_xlim((0.5,1.5))
#ax2.set_label('Blah')
#cbar.ax.text(16, 105., 'Phase Deviation (degree)', fontsize=18, va='center')
plt.savefig('HISTO_amp.png', format='PNG', dpi=400)


