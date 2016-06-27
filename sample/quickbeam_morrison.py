from netCDF4 import Dataset
import numpy as np
import quickbeam
import sys
import misc.fileops
from time import *
import os.path
from multiprocessing import Pool
from glob import glob

def create_outfile(inputfile):
    return '/home/cumulus/cp/Users/gryspeerdt/mor_'+os.path.basename(inputfile)+'.refl2m'


def calc_refl(inputfile):
    outputfile = create_outfile(inputfile)
    print os.path.basename(outputfile)
    #test = '/home/edward/LocalData/wrf_testdata.nc'

    qb = quickbeam.Quickbeam()
    qb.settings['hclass_file'] = 'data/hclass_wrf_morrison.dat'
    qb.read_hclass()

    ncdf = Dataset(inputfile)
    refl_att = np.zeros(ncdf.variables['TT'].shape)
    refl = np.zeros(ncdf.variables['TT'].shape)

    for time in range(refl.shape[0]):
        for lon in range(refl.shape[-1]):
            if lon%100 == 0:
                print time, lon, strftime("%a, %d %b %Y %H:%M:%S"), os.path.basename(outputfile)
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
            qb.hclass[1]['p1'] = ncdf.variables[
                'QNRAIN'][time, ::-1, lon, lat_min:lat_max]

            qb.hclass[2]['data'] = ncdf.variables[
                'QSNOW'][time, ::-1, lon, lat_min:lat_max] * 1000
            qb.hclass[2]['p1'] = ncdf.variables[
                'QNSNOW'][time, ::-1, lon, lat_min:lat_max]

            qb.hclass[3]['data'] = ncdf.variables[
                'QGRAUP'][time, ::-1, lon, lat_min:lat_max] * 1000
            qb.hclass[3]['p1'] = ncdf.variables[
                'QNGRAUPEL'][time, ::-1, lon, lat_min:lat_max] * 1000

            qb.hclass[4]['data'] = ncdf.variables[
                'QICE'][time, ::-1, lon, lat_min:lat_max] * 1000
            qb.hclass[4]['p1'] = ncdf.variables[
                'QNICE'][time, ::-1, lon, lat_min:lat_max]

            res = qb.radarsim()
            refl_att[time, ::-1, lon, lat_min:lat_max] = res['Zcorr']
            refl[time, ::-1, lon, lat_min:lat_max] = res['Zeff']

    ncdf.close()

    misc.fileops.nc4_dump_mv(outputfile, {'REFL_94GHz_att':refl_att,
                                          'REFL_94GHz':refl})


if __name__ == '__main__':
    pool = Pool(20)
    files = glob('/home/cumulus/cp/Users/white/data/wrf_congo_basin/ERA-Interim/archer_4km_august_2007_fullrun_pbl_scheme/pre-interp/heightlevs/interp/*INTRP')
    pool.map(calc_refl, files)
