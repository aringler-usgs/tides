#!/usr/bin/env python
from obspy.core import Stream, Trace, UTCDateTime
import os
from obspy.clients.fdsn import Client
import numpy as np
import glob
from scipy.signal import hilbert
from obspy.signal.cross_correlation import correlate, xcorr_max
import matplotlib.pyplot as plt


import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

def get_fb(tidetype):
    if tidetype == "all":
        fm = 1./(30.*60.*60.)
        fM = 1./(11.*60.*60.)
    elif tidetype == "semidirunal":
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


def get_syns_data(net, sta, loc, stime, etime, srate, tidetype, comps):
    client = Client('IRIS')
    inv = client.get_stations(network=net, station = sta, location=loc, starttime=stime,
                              endtime = etime, channel = comps, level="response")

    coord = inv.get_coordinates(net + '.' + sta + '.00.LHZ', stime)

    st = spotlme(coord['latitude'], coord['longitude'], coord['elevation'], False, stime, etime, srate)
    stCorr = spotlme(coord['latitude'], coord['longitude'], coord['elevation'], True, stime, etime, srate)
    for chan in ['Z', 'N', 'E']:
        trT = stCorr.select(component = chan)[0]
        #print(trT)
        #print(st.select(component = chan)[0])
        trT.data += st.select(component = chan)[0].data[:]
        st += trT
    # Remove all non-load corrected synthetics
    for tr in st.select(network="YY"):
        st.remove(tr)
    st.detrend('constant')
    # Remove fake response to avoid processing artifacts in the synthetics
    paz = {'zeros': [0. -1.j], 'poles': [0. -1.j], 'sensitivity': 1., 'gain': 1.}
    st.simulate(paz_remove=paz, taper=False)
    st_data = client.get_waveforms(net, sta, loc, comps, stime, etime)
    st_data.detrend('constant')
    st_data.merge(fill_value=0.)
    st_data.remove_response(inv, 'ACC', taper=False)
    for tr in st_data:
        tr.data *= 10**8
    if comps != 'LHZ':
        st_data.rotate('->ZNE', inventory=inv)
    st += st_data
    return st

def proc_tides(ctime, net, sta, loc):
    # Helper function for windowing and calculation
    stime = ctime -5*24*60*60
    etime = ctime +5*24*60*60
    srate = 1
    tidetype = 'semidirunal'
    comps ='LHZ'
    st = get_syns_data(net, sta, loc, stime, etime, srate, tidetype, comps)
    fm, fM = get_fb(tidetype)
    st.filter('bandpass',freqmin=fm, freqmax=fM, zerophase=True, corners=2)
    st.sort()
    # we now have everything and can do the calulation
    #st.trim(ctime-5*24*60*60, ctime + 6*24*60*60)
    st2 = st.select(component='Z')

    cc = correlate(st2[0].data, st2[1].data, 1000)
    shift, val = xcorr_max(cc)
    ptp = np.ptp(st2[0].data)
    ptp2 = np.ptp(st2[1].data)
    return shift, val, ptp, ptp2, st2

net, sta, loc = 'IU', 'ANMO', '00'

#f = open(net +'_' + sta +'_' + loc +'results','w')
#f.write('Phase Shift, Correlation, ptp-data, ptp-synthetic, DOY, Year\n')

ctime = UTCDateTime('2019-001T00:00:00')
etime = UTCDateTime('2020T001T00:00:00')

# while ctime <= etime:
#     try:
#         shift, val, ptp, ptp2, st2  = proc_tides(ctime, net, sta, loc)
#         print(str(shift) + ' ' + str(val))
#         #f.write(str(shift) +', ' + str(val) + ', ' + str(ptp) + ', ' + str(ptp2) + ', ' + 
#         #        str(ctime.julday) + ', ' + str(ctime.year) + '\n')
#     except:
#         print('problem')
#     ctime += 24*60*60

def mktime_series(D, psi, lon, lat, t):
    # eqn 1 of Ray et al. (2001)
    phi = (lat + 90.)*np.pi/180.
    lon *= np.pi/180.
    w = 2.*np.pi/(12.*60.*60. + 25.*60.)
    zeta = -D*np.cos(w*t + 2*phi - psi)*3*np.sin(lon)**2
    return zeta






shift, val, ptp, ptp2, st2  = proc_tides(ctime, net, sta, loc)







fig = plt.figure(1)
for idx, comp in enumerate(['Z']):
    for idx, tr in enumerate(st2):
        plt.plot(tr.times()/(24*60*60),tr.data, label=tr.id + ' ' + str(shift) + ' s')
        zeta = mktime_series(20, 0.,-106.4572, 34.94591,tr.times() )
        plt.plot(tr.times()/(24*60*60), zeta, label='Harmonic')
    plt.xlim((min(tr.times()/(24*60*60)),max(tr.times()/(24*60*60))))
    plt.ylabel('Acceleration ($\mu$Gal)')
    plt.legend(prop={'size': 12}, loc=2)
plt.xlabel('Time (days)')
plt.show()