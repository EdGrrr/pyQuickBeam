from netCDF4 import Dataset
import numpy as np
import quickbeam
import pkg_resources
import misc.fileops

inputfile = '/home/edward/LocalData/wrf_testdata.nc'

mphys_parametrisation = 'LIN' # Other choices are 'THOMPSON' or 'MORRISON'

qb = quickbeam.Quickbeam()
if mphys_parametrisation == 'LIN':
    qb.settings['hclass_file'] = pkg_resources.resource_filename(
        'pyQuickBeam', 'data/hclass_wrf.dat')
elif mphys_parametrisation == 'THOMPSON':
    qb.settings['hclass_file'] = pkg_resources.resource_filename(
        'pyQuickBeam', 'data/hclass_wrf_thompson.dat')
elif mphys_parametrisation == 'MORRISON':
    qb.settings['hclass_file'] = pkg_resources.resource_filename(
        'pyQuickBeam', 'data/hclass_wrf_morrison.dat')
else:
    print('Parametrisation not implemented')
    exit()
qb.read_hclass()

with Dataset(inputfile) as ncdf:
    refl_att = np.zeros(ncdf.variables['TT'].shape)
    refl = np.zeros(ncdf.variables['TT'].shape)

    for time in range(refl.shape[0]):
        for lon in range(refl.shape[-1]):
            if lon%100 == 0:
                print(time, lon)
            lat_min, lat_max = 0, refl.shape[-1]

            qb.met['Height'] = np.arange(
                20.5, 0, -1)[:, None].repeat(lat_max - lat_min, axis=1)
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
            if mphys_parametrisation in ['THOMPSON', 'MORRISON']:
                qb.hclass[1]['p1'] = ncdf.variables[
                    'QNRAIN'][time, ::-1, lon, lat_min:lat_max]

            qb.hclass[2]['data'] = ncdf.variables[
                'QSNOW'][time, ::-1, lon, lat_min:lat_max] * 1000
            if mphys_parametrisation == 'MORRISON':
                qb.hclass[2]['p1'] = ncdf.variables[
                    'QNSNOW'][time, ::-1, lon, lat_min:lat_max]

            qb.hclass[3]['data'] = ncdf.variables[
                'QGRAUP'][time, ::-1, lon, lat_min:lat_max] * 1000
            if mphys_parametrisation == 'MORRISON':
                qb.hclass[3]['p1'] = ncdf.variables[
                    'QNGRAUP'][time, ::-1, lon, lat_min:lat_max] * 1000
            elif mphys_parametrisation == 'THOMPSON':
                qb.hclass[3]['p1'] = ncdf.variables[
                    'QRAIN'][time, ::-1, lon, lat_min:lat_max] * 1000
                qb.hclass[3]['p2'] = ncdf.variables[
                    'QNRAIN'][time, ::-1, lon, lat_min:lat_max]

            qb.hclass[4]['data'] = ncdf.variables[
                'QICE'][time, ::-1, lon, lat_min:lat_max] * 1000
            if mphys_parametrisation in ['THOMPSON', 'MORRISON']:
                qb.hclass[4]['p1'] = ncdf.variables[
                    'QNICE'][time, ::-1, lon, lat_min:lat_max]

            res = qb.radarsim()
            refl_att[time, ::-1, lon, lat_min:lat_max] = res['Zcorr']
            refl[time, ::-1, lon, lat_min:lat_max] = res['Zeff']

misc.fileops.nc4_dump_mv(inputfile+'.refl', {'REFL_94GHz_att':refl_att,
                                             'REFL_94GHz':refl})
