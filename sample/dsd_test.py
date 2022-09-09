import pyQuickBeam
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
bpm = 2.
Q = 1.34

N1 = pyQuickBeam.dsd(Q,  # Mixing ratio
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
                     bpm)  # bpm

N2 = pyQuickBeam.dsd(Q,  # Mixing ratio
                     D,  # Size dist
                     1,  # Dist type
                     1.1,  # Density of air
                     -30,  # Temp
                     dmin,  # dmin
                     dmax,  # dmax
                     250000,
                     -2,
                     -3,
                     fc,  # fc
                     0,  # scaled?
                     apm,  # apm
                     bpm)  # bpm

N3 = pyQuickBeam.dsd(Q,  # Mixing ratio
                     D,  # Size dist
                     7,  # Dist type
                     1.1,  # Density of air
                     -30,  # Temp
                     dmin,  # dmin
                     dmax,  # dmax
                     Q*100,
                     10000,
                     -1,
                     fc,  # fc
                     0,  # scaled?
                     apm,  # apm
                     bpm)  # bpm

for dsdN, name in [(N1, 'Field 2005'),
                   (N2, 'Modified Gamma'),
                   (N3, 'Thompson Graupel')]:
    print('{: >20}. Num:{:.0f} Mass: {:.0f}'.format(
        name,
        (Dwidth*dsdN).sum()*1e6/1.071,
        (Dwidth*dsdN*apm*(Dmean)**bpm).sum()))

    plt.plot(D, dsdN, label=name)

plt.semilogx()
plt.semilogy()

plt.xlim(D[0], D[-1])
plt.ylim(1e-14, 1e-1)

plt.xlabel('D (um)')
plt.ylabel(r'$\frac{d\ln N}{d \ln D}$')

plt.legend()
plt.show()
