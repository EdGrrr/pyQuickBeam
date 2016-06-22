import radsim
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

# Constants
maxhclass = 50  # Number of hm classes
nd = 85  # Number of diameters
dmin, dmax = 0.1, 10000.  # Min and max particle size
D = np.exp(np.linspace(np.log(dmin), np.log(dmax), 85))
Dwidth = np.diff(D)
Dwidth = np.concatenate((Dwidth, np.array([Dwidth[-1]])))
Dmean = 0.5 * (D[1:]+D[:-1])
Dmean = np.concatenate((Dmean, np.array([Dmean[-1]])))
Dmean = D
fc = np.zeros(D.shape)

apm = (6/np.pi) * 100
bpm = 3.
Q = 1.34e-2

N1 = radsim.dsd(Q*100,  # Mixing ratio
                D,  # Size dist
                6,  # Dist type
                1.1,  # Density of air
                -30,  # Temp
                dmin,  # dmin
                dmax,  # dmax
                -1,
                -1,
                -1,
                fc,  # fc
                0,  # scaled?
                apm,  # apm
                2)  # bpm

N2 = radsim.dsd(Q,  # Mixing ratio
                D,  # Size dist
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

N3 = radsim.dsd(Q,  # Mixing ratio
                D,  # Size dist
                7,  # Dist type
                1.1,  # Density of air
                -10,  # Temp
                dmin,  # dmin
                dmax,  # dmax
                Q*100,
                10000,
                -1,
                fc,  # fc
                0,  # scaled?
                apm,  # apm
                bpm)  # bpm
print 'Field 2005 distn'
print 'N1 num : ', (Dwidth*N1).sum()*1e6/1.071
print 'N1 mass: ', (Dwidth*N1*apm*(Dmean)**2).sum()*1e-1
print ''
print 'Modified Gamma'
print 'N2 num : ', (Dwidth*N2).sum()*1e6/1.071
print 'N2 mass: ', (Dwidth*N2*apm*(Dmean)**bpm).sum()*1e-1
print ''
print 'Thompson graupel'
print 'N3 num : ', (Dwidth*N3).sum()*1e6/1.071
print 'N3 mass: ', (Dwidth*N3*apm*(Dmean)**bpm).sum()*1e-1


plt.plot(np.log(D), np.log(N1))
plt.plot(np.log(D), np.log(N2))
plt.plot(np.log(D), np.log(N3))
plt.show()
