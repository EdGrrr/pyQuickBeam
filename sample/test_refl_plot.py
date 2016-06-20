from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import quickbeam
import sys
import misc.fileops
from time import *

#inputfile = sys.argv[1]
#outputfile = sys.argv[2]

inputfile = '/home/edward/LocalData/wrf_testdata.nc'

qb = quickbeam.Quickbeam()
qb.settings['hclass_file'] = 'data/hclass_wrf.dat'
qb.read_hclass()

ncdf = Dataset(inputfile)
refl_att = np.zeros(ncdf.variables['TT'].shape)
refl = np.zeros(ncdf.variables['TT'].shape)

time,lon = 0,250
lat_min, lat_max = 0, refl.shape[-1]

#for i in [1,2,3,4]:
#    qb.hclass[i]['type'] = -1


qb.met['Height'] = np.arange(
    20.5, 0, -1)[:, None].repeat(lat_max - lat_min, axis=1) * 1000
qb.met['Pressure'] = ncdf.variables['PRES'][
    time, ::-1, lon, lat_min:lat_max] / 100.
qb.met['Temperature'] = ncdf.variables[
    'TT'][time, ::-1, lon, lat_min:lat_max]
qb.met['RH'] = ncdf.variables[
    'RH'][time, ::-1, lon, lat_min:lat_max]
qb.hclass[0]['data'] = ncdf.variables[
    'QCLOUD'][time, ::-1, lon, lat_min:lat_max] * 1000

qb.hclass[1]['data'] = ncdf.variables[
    'QRAIN'][time, ::-1, lon, lat_min:lat_max] * 1000
qb.hclass[1]['number'] = ncdf.variables[
    'QNRAIN'][time, ::-1, lon, lat_min:lat_max]

qb.hclass[2]['data'] = ncdf.variables[
    'QSNOW'][time, ::-1, lon, lat_min:lat_max] * 1000

qb.hclass[3]['data'] = ncdf.variables[
    'QGRAUP'][time, ::-1, lon, lat_min:lat_max] * 1000

qb.hclass[4]['data'] = ncdf.variables[
    'QICE'][time, ::-1, lon, lat_min:lat_max] * 1000
qb.hclass[4]['number'] = ncdf.variables[
    'QNICE'][time, ::-1, lon, lat_min:lat_max]

res = qb.radarsim()

#ncdf.close()

plt.subplot(311)
plt.imshow(res['Zeff'],vmin=-30,vmax=30)
plt.title('Refl two moment Rain and Ice')
plt.colorbar()

qb.hclass[0]['p3'] = 2
qb.hclass[1]['p1'] = -1
qb.hclass[1]['p2'] = 1000
qb.hclass[4]['p1'] = 100
qb.hclass[4]['p2'] = -1

res2 = qb.radarsim()

plt.subplot(312)
plt.imshow(res2['Zeff'],vmin=-30,vmax=30)
plt.title('Refl crappy single moment')
plt.colorbar()

plt.subplot(313)
plt.imshow(res2['Zeff']-res['Zeff'],vmin=-3,vmax=3)
plt.title('Difference')
plt.colorbar()


plt.show()


