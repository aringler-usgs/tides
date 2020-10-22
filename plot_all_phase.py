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
    lats, lons, amps_good, amps_goodstd, lens_good = {}, {}, {}, {}, {}

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
                    phase = float(result[6])*180/np.pi
                    if phase > 10:
                        phase = 10.
                    if phase < -10.:
                        phase = -10.
                    amps.append(phase)
            lens_good[sta + '_' + comp] = len(amps)
            amps_good[sta + '_' + comp] = np.mean(amps)
            amps_goodstd[sta + '_' + comp] = np.std(amps)
        try:
        #if True:
            inv = client.get_stations(starttime=stime, endtime=etime, network=net,
                 station = sta, channel="LHZ", level="response")
            if sta == 'GRFO':
                coors = inv.get_coordinates(net +'.' + sta + '..LHZ', stime)
            else:
                coors = inv.get_coordinates(net +'.' + sta + '.00.LHZ', stime)
        except:
            print('BAHAHAHAHAH' + sta)
            continue
        lats[sta] = coors['latitude']
        lons[sta] = coors['longitude']


    letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
    for idx, comp in enumerate(['Z', 'N', 'E']):
        vals, std, alats, alons, lens = [], [], [], [], []
        for sta in stas:
            try:
                vals.append(amps_good[sta + '_' + comp])
                std.append(amps_goodstd[sta + '_' + comp])
                alats.append(lats[sta])
                alons.append(lons[sta])
                lens.append(100*lens_good[sta + '_' + comp]/365.)
            except:
                continue

        ax = fig.add_subplot(3,2,2*idx+idx2+ 1, projection = ccrs.Robinson())
        ax.coastlines()
        ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
        ax.set_global()
        ax.coastlines()
        ax.set_title(letters[2*idx + idx2 ] + ' ' + comp +  ' ' + tidetype.replace('_',''), loc='left')
        im = ax.scatter(alons, alats, c=vals, transform=ccrs.Geodetic(), zorder=3, vmin=-10, vmax = 10)
        

cb_ax = fig.add_axes([0.2, 0.05, 0.6, 0.03])
cbar = fig.colorbar(im, cax=cb_ax, orientation='horizontal')
cbar.set_label('Phase Difference (degree)') 
#pos = cbar.ax.get_position()
#ax2 = cbar.ax.twiny()
#ax2.set_xlim((-10.,10.))
#ax2.set_xlim((0.5,1.5))
#ax2.set_label('Blah')
#cbar.ax.text(16, 105., 'Phase Deviation (degree)', fontsize=18, va='center')
plt.savefig('Map_PHASE_NEW.png', format='PNG', dpi=400)

# def setupmap(handle):


#     handle.add_feature(cfeature.LAND)
#     handle.add_feature(cfeature.OCEAN)
#     handle.add_feature(cfeature.COASTLINE)
#     handle.add_feature(cfeature.BORDERS, linestyle=':')
#     handle.add_feature(cfeature.LAKES, alpha=0.5)
#     handle.add_feature(cfeature.RIVERS)
#     handle.add_feature(cfeature.STATES, edgecolor='gray')
#     return handle
# letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
# fig = plt.figure(1, figsize=(12,14))
# for idx, f0 in enumerate(f0s):
#     ap,astd, alat, alon = [], [], [], []
#     ax = fig.add_subplot(3,2,idx+1, projection = ccrs.Robinson())
#     ax.set_title(letters[idx] + ' ' + str(int(round(1./f0,0))) + ' s', loc='left')
#     ax.coastlines()
#     ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
#     ax.set_global()
#     ax.coastlines()
#     #ax=setupmap(ax)
#     for ctrip in trip:
#         try:
#         #if True:
#             net, sta = ctrip.split('_')
#             camps = np.array(amps[sta])
#             evtemp = np.array(events[sta])
#             cf0sused = np.array(f0sused[sta])
#             camps = camps[(cf0sused == f0)]
#             camps = camps[~np.isnan(camps)]
#             evevec = evtemp[(cf0sused == f0)]
#             num = 100*sum(evevec)/len(evevec)
#             if len(camps) < 2:
#                 continue
#             ap.append(num)
#             #astd.append(np.std(camps))
#             alat.append(lats[sta])
#             alon.append(lons[sta])
#             print(sta + ' ' + str(num) + ' ' + str(f0))
            
#         except:
#             continue
    
    

#     im = ax.scatter(alon, alat, c=ap, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=100.)


# cb_ax = fig.add_axes([0.2, 0.05, 0.6, 0.03])
# cbar = fig.colorbar(im, cax=cb_ax, orientation='horizontal')

# cbar.set_label('Percentage Used (\%)') 

# 
# plt.show()
