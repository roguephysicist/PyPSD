'''
usage: python pypsd.py [-h] inputfile binsfile

A Python script for calculating the particle size distribution (PSD) of any
sample.

Written by Sean M. Anderson and Liliana Villafana-Lopez
'''

import argparse
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Argparse stuff
PARSER = argparse.ArgumentParser()
PARSER.add_argument("inputfile", help="File containing the particle areas", type=str)
PARSER.add_argument("binsfile", help="File containing the particle diameter bins", type=str)
ARGS = PARSER.parse_args()

# Read input file from command line, create arrays
INFILE = ARGS.inputfile
BASENAME = INFILE.split('.txt')[0]
PARTICLES = np.loadtxt(INFILE)
BINS = np.loadtxt(ARGS.binsfile, unpack=True)

def distribution(values, cutoff):
    """ a function for creating and then solving a linear equation """
    counter = np.argmax(values >= cutoff)
    point2 = np.array([BINS[counter], values[counter]])
    point1 = np.array([BINS[counter-1], values[counter-1]])
    slope = (point2[1] - point1[1])/(point2[0] - point1[0])
    intercept = point2[1] - slope*point2[0]
    dist = (cutoff - intercept) * (1/slope)
    return dist

# Rolling mean for significant sample size
AVGCUM = np.cumsum(PARTICLES, dtype=float)/np.arange(1, PARTICLES.size+1)

# Binning and frequencies
DIAMETERS = 2 * np.sqrt(PARTICLES/np.pi)
NI = np.bincount(np.digitize(DIAMETERS, BINS), minlength=BINS.size) # dmin
AI = 4 * np.pi * (BINS/2)**2 * NI
VI = (4.0/3.0) * np.pi * (BINS/2)**3 * NI

# Differential percentage
NPCT = NI * (1/np.sum(NI)) * 100
APCT = AI * (1/np.sum(AI)) * 100
VPCT = VI * (1/np.sum(VI)) * 100

# Cumulative percentage
NCUM = np.cumsum(NPCT, dtype=float)
ACUM = np.cumsum(APCT, dtype=float)
VCUM = np.cumsum(VPCT, dtype=float)

# Data for distributions
np.savetxt(BASENAME + '_distdata.txt', \
           np.c_[BINS, NPCT, APCT, VPCT], \
           fmt=('%07.3f', '%07.3f', '%07.3f', '%07.3f'), delimiter='    ', \
           header='Bins[um]'+1*" "+ 'Number'+5*" "+'Area'+7*" "+'Volume')

### Plot with matplotlib
FIG = plt.figure(1, figsize=(12, 9))
FIG.set_tight_layout(True)

# Upper subplot, significant samples
plt.subplot(211)
plt.grid(True, which="both")
plt.title("Significant Sample Average")
plt.xlabel("Samples")
plt.ylabel("Cumulative Average")
plt.plot(AVGCUM, lw=2.0)

# Lower subplot, size distributions
AX = plt.subplot(212)
plt.title("Droplet Size Distribution")
plt.xlabel("Diameter (um)")
plt.ylabel("Differential (%)")
plt.grid(True, which="both")
plt.xlim((0.01, 100.0))
plt.xscale('log')
AX.xaxis.set_major_formatter(ScalarFormatter())
plt.plot(BINS, NPCT, label='Number', lw=2.5, color='red')
plt.plot(BINS, APCT, label='Area', lw=2.5, color='purple')
plt.plot(BINS, VPCT, label='Volume', lw=2.5, color='green')
plt.legend()

# Save plot to .png file
FIG.savefig(BASENAME + '_distributions.png')

# Number distributions
ND10 = distribution(NCUM, 10)
ND50 = distribution(NCUM, 50)
ND90 = distribution(NCUM, 90)

# Area distributions
AD10 = distribution(ACUM, 10)
AD50 = distribution(ACUM, 50)
AD90 = distribution(ACUM, 90)

# Volume distributions
VD10 = distribution(VCUM, 10)
VD50 = distribution(VCUM, 50)
VD90 = distribution(VCUM, 90)

# Span
NSPAN = (ND90 - ND10)/ND50
ASPAN = (AD90 - AD10)/AD50
VSPAN = (VD90 - VD10)/VD50

# Mode
NMODE = BINS[np.argmax(NPCT)]
AMODE = BINS[np.argmax(APCT)]
VMODE = BINS[np.argmax(VPCT)]

# Median
NMEDIAN = BINS[np.argmax(NCUM >= 50)]
AMEDIAN = BINS[np.argmax(ACUM >= 50)]
VMEDIAN = BINS[np.argmax(VCUM >= 50)]

# D[1,0], D[3,2], D[4,3]
D_1_0 = np.sum(NI * BINS)/np.sum(NI)
D_3_2 = np.sum(NI * BINS**3)/np.sum(NI * BINS**2)
D_4_3 = np.sum(NI * BINS**4)/np.sum(NI * BINS**3)

# Print results to file
with open(BASENAME + '_granulometry.txt', 'w') as outfile:
    print(9*" " + 'Number' + 2*" " + 'Area' + 4*" " + 'Volume', file=outfile)
    print('D10:    {0:6.3f}  {1:6.3f}  {2:6.3f}'.format(ND10, AD10, VD10), file=outfile)
    print('D50:    {0:6.3f}  {1:6.3f}  {2:6.3f}'.format(ND50, AD50, VD50), file=outfile)
    print('D90:    {0:6.3f}  {1:6.3f}  {2:6.3f}'.format(ND90, AD90, VD90), file=outfile)
    print('Span:   {0:6.3f}  {1:6.3f}  {2:6.3f}'.format(NSPAN, ASPAN, VSPAN), file=outfile)
    print('Mode:   {0:6.3f}  {1:6.3f}  {2:6.3f}'.format(NMODE, AMODE, VMODE), file=outfile)
    print('Median: {0:6.3f}  {1:6.3f}  {2:6.3f}\n'.format(NMEDIAN, AMEDIAN, VMEDIAN), file=outfile)
    print('D[1,0]: {0:6.3f}'.format(D_1_0), file=outfile)
    print('D[3,2]: {0:6.3f}'.format(D_3_2), file=outfile)
    print('D[4,3]: {0:6.3f}\n'.format(D_4_3), file=outfile)
    print('Total particles: {0}'.format(PARTICLES.size), file=outfile)
