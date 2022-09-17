import matplotlib.pyplot as plt
import pyQuickBeam

qb = pyQuickBeam.Quickbeam()
qb.settings['do_ray'] = 1
refl = qb.radarsim()

plt.plot(refl['Zeff'], qb.met['Height'], label='Effective')
plt.plot(refl['Zcorr'], qb.met['Height'], label='Corrected')
plt.xlim(-20,10)
plt.ylim(0,10)
plt.xlabel('Reflectivity (dBZ)')
plt.ylabel('Height (km)')
plt.legend()
plt.show()

# f = open('../tests/radar_sim.txt', 'w')
# f.write('Height,Zeff,Zray,Zcorr,,\n')
# for i in range(len(qb.met['Height'])):
#     f.write(
#         '{},{},{},{},,\n'.format(
#             qb.met['Height'][i][0],
#             refl['Zeff'][i][0],
#             refl['Zray'][i][0],
#             refl['Zcorr'][i][0]))
# f.close()
