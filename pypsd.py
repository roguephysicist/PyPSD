'''
usage: python pypsd.py [-h] inputfile binsfile

A Python script for calculating the particle size distribution (PSD) of any
sample.

Written by Sean Anderson (https://github.com/roguephysicist/PyPSD)
'''

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter
import argparse

# Argparse stuff
parser = argparse.ArgumentParser()
parser.add_argument("inputfile", help="File containing the particle areas", type=str)
parser.add_argument("binsfile", help="File containing the diameter bins", type=str)
args = parser.parse_args()

# Read input file from command line, create arrays
infile = args.inputfile
basename = infile.split('.')[0]
particles = np.loadtxt(infile)
dmin, dmax, bins = np.loadtxt(args.binsfile, unpack=True)

# A function for creating and then solving a linear equation
def distribution(values, cutoff):
    n = np.argmax(values>=cutoff)
    point2 = np.array([bins[n], values[n]])
    point1 = np.array([bins[n-1], values[n-1]])
    slope = (point2[1] - point1[1])/(point2[0] - point1[0])
    intercept = point2[1] - slope*point2[0]
    dist = (cutoff - intercept) * (1/slope)
    return dist

# Rolling mean for significant sample size
avgcum = np.cumsum(particles, dtype=float)/np.arange(1, particles.size+1)

# Binning and frequencies
diameters = 2 * np.sqrt(particles/np.pi)
ni = np.bincount(np.digitize(diameters, bins), minlength=bins.size) # dmin
ai = 4 * np.pi * (bins/2)**2 * ni
vi = (4.0/3.0) * np.pi * (bins/2)**3 * ni

# Differential percentage
npct = ni * (1/np.sum(ni)) * 100
apct = ai * (1/np.sum(ai)) * 100
vpct = vi * (1/np.sum(vi)) * 100

# Cumulative percentage
ncum = np.cumsum(npct, dtype=float)
acum = np.cumsum(apct, dtype=float)
vcum = np.cumsum(vpct, dtype=float)

# Data for distributions
np.savetxt(basename + '_distdata.txt', \
           np.c_[bins, npct, apct, vpct], \
           fmt=('%07.3f', '%07.3f', '%07.3f', '%07.3f'), delimiter='    ', \
           header='Bins[um]'+1*" "+ 'Number'+5*" "+'Area'+7*" "+'Volume')

### Plot with matplotlib
fig = plt.figure(1, figsize=(12, 9))
fig.set_tight_layout(True)

# Upper subplot, significant samples
plt.subplot(211)
plt.grid(True, which="both")
plt.title("Significant Sample Average")
plt.xlabel("Samples")
plt.ylabel("Cumulative Average")
plt.plot(avgcum, lw=2.0)

# Lower subplot, size distributions
ax = plt.subplot(212)
plt.title("Droplet Size Distribution")
plt.xlabel("Diameter (um)")
plt.ylabel("Differential (%)")
plt.grid(True, which="both")
plt.xlim((0.01,100.0))
plt.xscale('log')
ax.xaxis.set_major_formatter(ScalarFormatter())
plt.plot(bins, npct, label='Number', lw=2.5, color='red')
plt.plot(bins, apct, label='Area', lw=2.5, color='purple')
plt.plot(bins, vpct, label='Volume', lw=2.5, color='green')
plt.legend()

fig.savefig(basename + '_distributions.png') # saves to png file
# plt.ion()
# plt.pause(0.001)
# plt.show()
###

# User input for peaks
# print('Are there any non-monomodal distributions? Just leave blank if not.')
# distro = input('(Options: number, area, volume): ') or 'none'
# if distro is not 'none':
#     print('\n{0} distribution selected.'.format(distro))
#     user_diams = input('Input approx values for mins, separated by commas: ')
#     a = [float(x) for x in user_diams.split(',')]
#     print(a)


# Number distributions
nD10 = distribution(ncum, 10)
nD50 = distribution(ncum, 50)
nD90 = distribution(ncum, 90)

# Area distributions
aD10 = distribution(acum, 10)
aD50 = distribution(acum, 50)
aD90 = distribution(acum, 90)

# Volume distributions
vD10 = distribution(vcum, 10)
vD50 = distribution(vcum, 50)
vD90 = distribution(vcum, 90)

# Span
nspan = (nD90 - nD10)/nD50 
aspan = (aD90 - aD10)/aD50 
vspan = (vD90 - vD10)/vD50 

# Mode
nmode = bins[np.argmax(npct)]
amode = bins[np.argmax(apct)]
vmode = bins[np.argmax(vpct)]

# Median
nmedian = bins[np.argmax(ncum>=50)]
amedian = bins[np.argmax(acum>=50)]
vmedian = bins[np.argmax(vcum>=50)]

# D[1,0], D[3,2], D[4,3]
D_1_0 = np.sum(ni * bins)/np.sum(ni)
D_3_2 = np.sum(ni * bins**3)/np.sum(ni * bins**2)
D_4_3 = np.sum(ni * bins**4)/np.sum(ni * bins**3)
# According to Dr. Villafana
# D_1_0 = np.sum(diameters)/np.sum(ni)
# D_3_2 = np.sum(diameters**3)/np.sum(diameters**2)
# D_4_3 = np.sum(diameters**4)/np.sum(diameters**3)

# Print results to file
with open(basename + '_granulometry.txt', 'w') as outfile:
    print(9*" " + 'Number' + 4*" " + 'Area' + 6*" " + 'Volume', file=outfile)
    print('D10:    {0:6.3f}    {1:6.3f}    {2:6.3f}'.format(nD10, aD10, vD10), file=outfile)
    print('D50:    {0:6.3f}    {1:6.3f}    {2:6.3f}'.format(nD50, aD50, vD50), file=outfile)
    print('D90:    {0:6.3f}    {1:6.3f}    {2:6.3f}'.format(nD90, aD90, vD90), file=outfile)
    print('Span:   {0:6.3f}    {1:6.3f}    {2:6.3f}'.format(nspan, aspan, vspan), file=outfile)
    print('Mode:   {0:6.3f}    {1:6.3f}    {2:6.3f}'.format(nmode, amode, vmode), file=outfile)
    print('Median: {0:6.3f}    {1:6.3f}    {2:6.3f}\n'.format(nmedian, amedian, vmedian), file=outfile)
    print('D[1,0]: {0:6.3f}'.format(D_1_0), file=outfile)
    print('D[3,2]: {0:6.3f}'.format(D_3_2), file=outfile)
    print('D[4,3]: {0:6.3f}\n'.format(D_4_3), file=outfile)
    print('Particles: {0}'.format(particles.size), file=outfile)

# plt.show()
