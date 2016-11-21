#!/Users/sma/anaconda3/bin/python
import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.ticker import ScalarFormatter
import sys

infile = sys.argv[1]
particles = np.loadtxt(infile)
dmin, dmax, bins = np.loadtxt("bins.dat", unpack=True)

def distribution(values, cutoff):
    n = np.argmax(values>cutoff)
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
# D_1_0 = np.sum(diameters)/np.sum(ni)
# D_3_2 = np.sum(diameters**3)/np.sum(diameters**2)
# D_4_3 = np.sum(diameters**4)/np.sum(diameters**3)

print('\nTotal number of particles: {0}\n'.format(particles.size))
print(9*" " + 'Number' + 4*" " + 'Surface' + 3*" " + 'Volume')
print(36*"=")
print('D10:    {0:6.3f}    {1:6.3f}    {2:6.3f}'.format(nD10, aD10, vD10))
print('D50:    {0:6.3f}    {1:6.3f}    {2:6.3f}'.format(nD50, aD50, vD50))
print('D90:    {0:6.3f}    {1:6.3f}    {2:6.3f}'.format(nD90, aD90, vD90))
print(36*"-")
print('Span:   {0:6.3f}    {1:6.3f}    {2:6.3f}'.format(nspan, aspan, vspan))
print('Mode:   {0:6.3f}    {1:6.3f}    {2:6.3f}'.format(nmode, amode, vmode))
print('Median: {0:6.3f}    {1:6.3f}    {2:6.3f}'.format(nmedian, amedian, vmedian))
print(36*"-")
print('D[1,0]: {0:6.3f}'.format(D_1_0))
print('D[3,2]: {0:6.3f}'.format(D_3_2))
print('D[4,3]: {0:6.3f}'.format(D_4_3))
print()

# plt.figure(1)

# plt.subplot(211)
# plt.grid(True)
# plt.title("Sample Average")
# plt.xlabel("Samples")
# plt.ylabel("Cumulative Average")
# plt.plot(avgcum, lw=2.0)

# plt.subplot(212)
# # ax.xaxis.set_major_formatter(ScalarFormatter())
# plt.title("Droplet Size Distribution")
# plt.xlabel("Diameter (um)")
# plt.ylabel("Differential (%)")
# plt.grid(True)
# # plt.xlim((0,10))
# plt.xscale('log')
# plt.plot(bins, npct, 'red', bins, apct, 'purple', bins, vpct, 'green', lw=2.5)

# # plt.tight_layout()
# plt.show()
