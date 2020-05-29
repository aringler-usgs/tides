#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)


sta = 'ANMO'
f = open('IU_' + sta + '_00results','r')
f.readline()
amps , ampsS =[],[]
for line in f:
    line = line.split(',')
    if float(line[1]) > 0.99:
        amps.append(float(line[2]))
        ampsS.append(float(line[3]))
f.close()

amps = np.array(amps)
ampsS = np.array(ampsS)
fig = plt.figure(1, figsize=(8,8))
plt.subplot(2,1,1)
plt.plot(amps, ampsS,'.')
plt.xlabel(sta + ' ($\mu$Gal)')
plt.ylabel('SPOTL ($\mu$Gal)')
plt.subplot(2,1,2)
plt.plot(amps/ampsS,'.')
plt.ylabel('Ratio')
plt.xlabel('Day of Year (2019)')
plt.savefig(sta + '.png', format='PNG')