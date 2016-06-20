import radsim
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

# Constants
maxhclass = 50  # Number of hm classes
nd = 85  # Number of diameters
dmin, dmax = 0.1, 10000.  # Min and max particle size
D = np.exp(np.linspace(np.log(dmin), np.log(dmax), 85))
fc = np.zeros(D.shape)

apm = (6/np.pi) * 100
bpm = 3
Q = 1.34e-2

N1 = radsim.dsd(Q,  # Mixing ratio
                D,  # Size dist
                # 100000,  # Number conc (-1 to use p values)
                1,  # Dist type
                1.1,  # Density of air
                -10,  # Temp
                dmin,  # dmin
                dmax,  # dmax
                250000,
                -2,
                -2,
                fc,  # fc
                0,  # scaled?
                apm,  # apm
                bpm)  # bpm
print (np.diff(D)*N1[:-1]).sum()*1e6/1.071
print (np.diff(D)*N1[:-1]*(0.5*(D[1:]+D[:-1]))**3).sum()*1e6

N2 = radsim.dsd(Q,  # Mixing ratio
                D,  # Size dist
                # 150000,  # Number conc (-1 to use p values)
                1,  # Dist type
                1.1,  # Density of air
                -10,  # Temp
                dmin,  # dmin
                dmax,  # dmax
                250000,
                -2,
                -3,
                fc,  # fc
                0,  # scaled?
                apm,  # apm
                bpm)  # bpm
print (np.diff(D)*N2[:-1]).sum()*1e6/1.071
print (np.diff(D)*N2[:-1]*(0.5*(D[1:]+D[:-1]))**3).sum()*1e6

D0 = 1e6*((Q*0.001*gamma(1))/
          (apm*250000*gamma(1+bpm)))**(1./bpm)
N3 = (1*250000/gamma(1))*(1e6/D0)*np.exp(-1*D/D0)*1e-12
print (np.diff(D)*N3[:-1]).sum()*1e6/1.071
print (np.diff(D)*N3[:-1]*(0.5*(D[1:]+D[:-1]))**3).sum()*1e6


plt.plot(np.log(D), np.log(N1))
plt.plot(np.log(D), np.log(N2))
plt.plot(np.log(D), np.log(N3))
plt.show()
