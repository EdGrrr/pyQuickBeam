import matplotlib.pyplot as plt
import quickbeam

qb = quickbeam.Quickbeam()
refl = qb.radarsim()

plt.plot(refl['Zeff'], qb.met['Height'], label='Effective')
plt.plot(refl['Zcorr'], qb.met['Height'], label='Corrected')
plt.xlim(-20,10)
plt.ylim(0,10)
plt.xlabel('Reflectivity (dBZ)')
plt.ylabel('Height (km)')
plt.legend()
plt.show()
