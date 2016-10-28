#!/Users/sma/anaconda3/bin/python
import numpy as np
import sys

infile = sys.argv[1]
particles = np.loadtxt(infile)
dmin, dmax, bins = np.loadtxt("bins.dat", unpack=True)

def cumulative(freq):
    pct = freq * (1/np.sum(freq)) * 100
    cum = np.cumsum(pct, dtype=float)
    return cum

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

ncum = cumulative(ni)
acum = cumulative(ai)
vcum = cumulative(vi)

print(6*" " + 'Number' + 10*" " + 'Surface' + 9*" " + 'Volume')
print(54*"-")
print('D10: {0:15.12f} {1:15.12f} {2:15.12f}'.format(distribution(ncum, 10), distribution(acum, 10), distribution(vcum, 10)))
print('D50: {0:15.12f} {1:15.12f} {2:15.12f}'.format(distribution(ncum, 50), distribution(acum, 50), distribution(vcum, 50)))
print('D90: {0:15.12f} {1:15.12f} {2:15.12f}'.format(distribution(ncum, 90), distribution(acum, 90), distribution(vcum, 90)))


# data = np.column_stack([bins, vcum])
# np.savetxt("out.dat", data)
