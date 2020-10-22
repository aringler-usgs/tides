#!/usr/bin/env python
from obspy.core import Stream, Trace, UTCDateTime
import os
from obspy.clients.fdsn import Client
import numpy as np
import os
import glob
from scipy.signal import hilbert
from obspy.signal.cross_correlation import correlate, xcorr_max
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.mlab import csd
#import subprocess
#subprocess.check_call(["latex"])
from obspy.signal.spectral_estimation import get_nhnm, get_nlnm

mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

client = Client('IRIS')

def get_fb(tidetype):
    if tidetype == "all":
        fm = 1./(30.*60.*60.)
        fM = 1./(11.*60.*60.)
    elif tidetype == "semidiurnal":
        fm = 1./(13.*60.*60.)
        fM = 1./(11.*60.*60.)
    elif tidetype == "diurnal":
        fm = 1./(30.*60.*60.)
        fM = 1./(22.*60.*60.)
    elif tidetype == "noise":
        fm = 1./(8.*60.*60.)
        fM = 1./(6.*60.*60.)
    elif tidetype == "dave":
        fm = 1./(4.*60.*60.)
        fM = 1./(3.*60.*60.)
    return fm, fM

def spotlme(lat, lon, dep, loading, stime, etime, srate, debug=True):
    if loading:

        tides = ['k1', 'k2', 'm2', 'm4', 'mf', 'mm', 'n2', 'o1', 'p1', 'q1', 's2']

        for idx, tide in enumerate(tides):
            string = '../bin/nloadf TEMP ' + str(lat) + ' ' + str(lon) + ' ' + str(dep) + ' '
            string += tide + '.osu.tpxo72.2010 ' + 'green.gbavap.std l'
            if idx == 0:
                pipe = ' >'
            else: 
                pipe = ' >>'
            cmd = string + ' ' + pipe + ' LoadALL' 
            os.system(cmd)
            if debug:
                print(cmd)
        st = Stream()
        for comp in ['Z', 'N', 'E']:
            if debug:
                print('On component:' + comp)
            
            cmd = 'cat LoadALL | ../bin/harprp '
            if comp == 'Z':
                cmd += 'g '
            elif comp == 'N':
                cmd += 't 0 '
            else:
                cmd += 't 90 '
            cmd += '> tempLoad2'    
            if debug:
                print(cmd)
            os.system(cmd)

            ctime = stime
            #cmd = 'cat tempLoad2 | ../bin/loadcomb t'
            tstring = str(ctime.year) + ' ' + str(ctime.julday) + ' 0 0 0'
            vals = int((etime -stime)/(srate))
            cmd = 'cat tempLoad2 | ../bin/hartid ' + tstring + ' ' + str(vals) + ' ' + str(srate) + ' > temp' + comp
            if debug:
                print(cmd)
            os.system(cmd)

            f=open('temp' + comp,'r')
            data = []
            for line in f:
                if comp == 'Z':
                    data.append(-float(line))
                else:
                    data.append(-float(line))

            f.close()
            data = np.array(data)
            stats = {'network': 'XX', 'station': sta, 'location': '',
                'channel' : 'LH' + comp, 'npts': len(data), 'sampling_rate': 1./srate,
                'mseed' : {'dataquality': 'D'}}
            stats['starttime'] = ctime
            st += Stream([Trace(data=data, header = stats)])
            ctime += 24.*60.*60.
            os.remove('temp' + comp)
            os.remove('tempLoad2')
        oldfiles = glob.glob('tempLoad.*')
        for curfile in oldfiles:
            os.remove(curfile)
        st.merge()
    else:
        # Here we just do the tides with no loading
        pfile = open('para_file','w')
        pfile.write(str(stime.year) + ',' + str(stime.julday) + ',0\n')
        pfile.write(str(etime.year) + ',' + str(etime.julday) + ',0\n')
        # In terms of 1 hour
        pfile.write(str(srate/(60*60)) + '\n')
        pfile.write('t\n')
        pfile.write(str(lat) + '\n')
        pfile.write(str(lon) + '\n')
        # compute gravity
        pfile.write('1\n')
        # Compute two tilt tides
        pfile.write('2\n')
        # compute no strain
        pfile.write('0\n')
        pfile.write('0\n')
        pfile.write('90\n')
        pfile.write('tempZ\n')
        pfile.write('tempN\n')
        pfile.write('tempE\n')
        pfile.close()
        cmd =  '../bin/ertid < para_file'
        os.system(cmd)
        os.remove('para_file')
        st = Stream()
        for comp in ['Z', 'N', 'E']:
            f=open('temp' + comp,'r')
            data = []
            for line in f:
                if comp == 'Z':
                    data.append(-float(line))
                else:
                    data.append(-float(line))

            f.close()
            os.remove('temp' + comp)
            data = np.asarray(data)
            stats = {'network': 'YY', 'station': sta, 'location': '',
                    'channel' : 'UH' + comp, 'npts': len(data), 'sampling_rate': 1./srate,
                    'mseed' : {'dataquality': 'D'}}
            stats['starttime'] = stime 
            st += Stream([Trace(data=data, header = stats)])
    
    return st


def get_syns_data(net, sta, loc, stime, etime, srate, tidetype):
    
    inv = client.get_stations(network=net, station = sta, location=loc, starttime=stime,
                              endtime = etime, channel = "LH*", level="response")
    if loc == '*':
        coord = inv.get_coordinates(net + '.' + sta + '.00.LHZ', stime)
    else:
        coord = inv.get_coordinates(net + '.' + sta + '.' + loc + '.LHZ', stime)

    st = spotlme(coord['latitude'], coord['longitude'], coord['elevation'], False, stime, etime, srate)
    stCorr = spotlme(coord['latitude'], coord['longitude'], coord['elevation'], True, stime, etime, srate)
    for chan in ['Z', 'N', 'E']:
        trT = stCorr.select(component = chan)[0]
        trHar = st.select(component = chan)[0]
        if trT.stats.npts > trHar.stats.npts:
            trT.data = trT.data[:trHat.stats.npts]
        if trT.stats.npts < trHar.stats.npts:
            trHar.data = trHar.data[:trT.stats.npts]

        #print(st.select(component = chan)[0])

        trT.data += trHar
        st += trT
    # Remove all non-load corrected synthetics
    for tr in st.select(network="YY"):
        st.remove(tr)
    st.detrend('constant')

    # Remove fake response to avoid processing artifacts in the synthetics
    paz = {'zeros': [0. -1.j], 'poles': [0. -1.j], 'sensitivity': 1., 'gain': 1.}
    st.simulate(paz_remove=paz, taper=False)
    st_data = client.get_waveforms(net, sta, loc, "LH*", stime, etime)
    st_data.detrend('constant')
    st_data.merge(fill_value=0.)
    st_data.remove_response(inv, 'ACC', taper=False)
    for tr in st_data:
        print(tr)
        tr.data *= 10**8

    st_data.rotate('->ZNE', inventory=inv)
    st += st_data
    return st

def cp(tr1, tr2, lenfft, lenol, delta):
    # Cross-power function
    sr = 1./float(delta)
    cpval,fre = csd(tr1.data, tr2.data, NFFT=lenfft, Fs=sr, noverlap=int(lenol*lenfft), scale_by_freq=True)
    fre = fre[1:]
    cpval = cpval[1:]
    return cpval, fre





net, loc= 'IU', '00'
stime = UTCDateTime('2019-050T00:00:00')
etime = UTCDateTime('2019-060T00:00:00')

inv = client.get_stations(network=net, station = '*', location=loc, starttime=stime,
                              endtime = etime, channel = "LH*")

stas = [sta.code for sta in inv[0]]

per_nlnm, pow_nlnm = get_nlnm()

per_nhnm, pow_nhnm = get_nhnm()
plt.semilogx(per_nhnm,pow_nhnm, linewidth=2, color='k', label='NLNM/NHNM')


stas = ['QSPA']
for sta in stas:
    tidetype = 'all'
    srate = 1

    
    st = get_syns_data(net, sta, loc, stime, etime, srate, tidetype)
    fm, fM = get_fb(tidetype)





    #st.filter('bandpass',freqmin=fm, freqmax=fM, zerophase=True, corners=2)
    st.sort()
    fig, axs = plt.subplots(3,1, figsize=(12,12))
    ax = axs.flatten()
    letts = ['(a)', '(b)', '(c)' ]
    for idx, comp in enumerate(['Z', 'N', 'E']):
        st2 = st.select(component=comp)
        if idx == 0:
            ax[0].set_title(sta + ' Start Time ' + str(stime.year) + ' DOY: ' + 
                str(stime.julday).zfill(3) + ' ' + str(stime.hour).zfill(2) +':' + 
                str(stime.minute).zfill(2) )
        for tr in st2:
            print(tr)
            p, f = cp(tr.data/(10**8),tr.data/(10**8), tr.stats.npts, 0, 1)


            if tr.stats.network == 'XX':
                ax[idx].semilogx((1/f)/(60*60),10*np.log10(p), label='Synthetic Tide Load Corrected')
                #ax[idx].plot(tr.times()/(24*60*60), tr.data, label='Synthetic Tide Load Corrected')
            elif tr.stats.network == 'YY':
                ax[idx].semilogx((1/f)/(60*60),10*np.log10(p), label='Synthetic Tide')
            else:
                ax[idx].semilogx((1/f)/(60*60),10*np.log10(p), label=(tr.id).replace('.',' '))
        ax[idx].set_xlim((9,33))
        ax[idx].text(8, -15, letts[idx])
        ax[idx].set_ylim((-140.,-20))
        #ax[idx].text(-1, 120., letts[idx])
        
        if idx == 2:
            ax[2].set_xlabel('Period (hour)')
        if idx == 1:
            ax[1].set_ylabel('PSD (dB rel. 1 $(m/s^2)^2/Hz$)')
        ax[idx].semilogx(per_nlnm/(60*60),pow_nlnm, linewidth=2, color='k')
        ax[idx].semilogx(per_nhnm/(60*60),pow_nhnm, linewidth=2, color='k', label='NLNM/NHNM')
        ax[idx].legend(ncol=3, loc=9, fontsize = 14)
    plt.savefig('PSD_PLOT.png')
                


     
       