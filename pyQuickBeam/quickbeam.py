import numpy as np
import struct
import pkg_resources
import radsim

DEBUG = False

metnames = ['Height', 'Pressure',
            'Temperature', 'RH']
hclasscols = ('id', 'type', 'col', 'phase',
              'cp', 'dmin', 'dmax',
              'apm', 'bpm', 'rho', 'p1', 'p2', 'p3', 'name')
settings = [('freq', 94),
            ('surface_radar', 0),
            ('use_mie_table', 0),
            ('use_gas_abs', 0),
            ('sonde_format', 2),
            ('do_ray', 0),
            ('melt_lay', 0),
            ('output_format', 2),
            ('output_disp', 0),
            ('k2', -1),
            ('hydro_file', pkg_resources.resource_filename(
                'pyQB', 'data/hydromet_example.dat')),
            ('hclass_file', pkg_resources.resource_filename(
                'pyQB', 'data/hclass.dat')),
            ('mie_table_name', pkg_resources.resource_filename(
                'pyQB', 'data/mie_table.dat'))]

# Constants
maxhclass = 50  # Number of hm classes
nd = 85  # Number of diameters
dmin, dmax = 0.1, 10000.  # Min and max particle size
D = np.exp(np.linspace(np.log(dmin), np.log(dmax), 85))


class Quickbeam:
    def __init__(self):
        '''Initialise a QuickBeam object. 

The object is populated with some defaults from the example files, so
you can run a simulation immediately with the 'radarsim' method.

The QuickBeam settings are:
  'freq': radar frequency
  'surface_radar': surface (1) or spaceborne (0) radar
  'use_mie_table': use pre-computed mie data (1)
  'use_gas_abs': calculate gas absorption (1)
  'do_ray': Perform rayleigh calculations (1)
  'k2': k2 value as default (-1)
  'hydro_file': hydrometeor data file (if using text file)
                load into Quickbeam instance with read_hydrodata()
                make sure to call read_class first if updating
                hydrometor properties
  'hclass_file': hydrometeor class data
                 load into Quickbeam instance with read_hclass()
  'mie_table_name': location of the pre-computed mie data.  Changing this
                    is not recommened as some assumptions on the table
                    structure are made in the code

Meteorological data is provided through the 'met' property. This data is
two dimensional and in the shape (ngate, nprofiles). The met data required
is:

  Height (m)
  Pressure (hPa)
  Temperature (K)
  RH (%)

Finally, hydrometeor data is provided through the 'hclass'
property. Each hclass gives the details of a hydrometeor type (size
distribution, mass-density relationship), along with the data required
to computer there in each gate. For each hydrometeor type, the
dictionary in hclass is

  'id': id number
  'type': size distribution function (-1 to ignore this class)
  'col': 
  'phase': Liquid (0) or Ice (1)
  'cp': TOCHECK - in quickbeam this is used to convert mixing ratios
        we are not doing that here, probably because g/kg assumed
  'dmin', 'dmax': Max and min of the size distribution,
                  currently ignored
  'apm', 'bpm': For calculating the density-size relationship
  'rho': Density (if using constant density
  'p1', 'p2', 'p3': Parameters for size distribution
  'name': Name of class (not required)
  'data': Hydrometeor concentration (g/ m3) 

The size distribution parameters can be cosntants, or arrays of model
data - for example, number concentration in 'p1' if using a modified
gamma (hclass type=1). 
        '''
        # The settings for running the code
        self.settings = {n[0]:n[1] for n in settings}
        # Now the hydrometeor data (put the data and parameters into dicts)
        self.read_hclass()
        self.read_hydrodata()
        # Read the Mie table
        self.read_mie_table()

    def radarsim(self):
        # Assume that whoever provided the hydrodata gave a sensible amount!
        # Input data is in the shape (ngate, nprofiles)
        rho_a = ((self.met['Pressure']*100) /
                 (287*(self.met['Temperature'])))
        ngate, nprofiles = self.met['Temperature'].shape
        z_vol = np.zeros(self.met['Temperature'].shape)
        z_ray = np.zeros(self.met['Temperature'].shape)
        kr_vol = np.zeros(self.met['Temperature'].shape)

        # Calculate frequency index for mie table
        if self.settings['use_mie_table'] == 1:
            freq_ind = np.argmin(abs(self.mt['freq']-self.settings['freq']))
            if min(abs(self.mt['freq']-self.settings['freq'])) > 3:
                raise(ValueError, 'Chosen frequency not in mie table')

        # For each hydrometeor
        for hm in self.hclass:
            if hm['type'] < 0:
                continue
            # Rho_eff calculation should only be done once
            #if (hm['phase'] == 1): now done for all types
            rho_eff = (6/np.pi)*hm['apm']*(D*1e-6)**(hm['bpm']-3)
            rho_eff[rho_eff < 5] = 5.
            rho_eff[rho_eff > 917] = 917.
            hm['f'] = np.digitize(rho_eff/917., self.mt['f_bins'])

            #Set twomoment flag
            tm_flag = (hm['p1'].shape == hm['data'].shape)
            if DEBUG:
                print(hm['name'], tm_flag)
            
            for g in range(z_vol.shape[0]):
                for p in range(z_vol.shape[1]):
                    if hm['data'][g, p] == 0:
                        continue
                    if self.settings['use_mie_table'] == 1:
                        # Identify the location in the met table
                        if hm['phase'] == 0:
                            temp_ind = int(
                                (self.met['Temperature'][g, p]-210.65)/5.)
                            temp_ind = np.clip(
                                temp_ind, 0, self.mt['cnt_liq']-1)
                            qext = self.mt['qext'][freq_ind, 0, temp_ind]
                            qbsca = self.mt['qbsca'][freq_ind, 0, temp_ind]
                        else:
                            temp_ind = int(
                                (self.met['Temperature'][g, p]-180.65)/5.)
                            temp_ind = np.clip(
                                temp_ind, 0, self.mt['cnt_ice']-1)
                            temp_ind += self.mt['cnt_liq']
                            qext = self.mt['qext'][freq_ind,
                                                   hm['f'],
                                                   temp_ind,
                                                   range(nd)]
                            qbsca = self.mt['qbsca'][freq_ind,
                                                     hm['f'],
                                                     temp_ind,
                                                     range(nd)]
                    else:
                        # Calculate the properties directly
                        qext, qbsca = radsim.optical_sphere(
                            self.settings['freq'],
                            D,
                            self.met['Temperature'][g, p]-273.15,
                            hm['phase'],
                            rho_eff)

                    if 'fc' not in hm.keys():
                        # DSD scaling factor - will be set by radsim.dsd
                        # if required
                        hm['fc'] = np.zeros(D.shape)
                        hm['scaled'] = 0

                    if hm['p1'].shape == hm['data'].shape:
                        p1 = hm['p1'][g, p]
                    else:
                        p1 = hm['p1']

                    if hm['p2'].shape == hm['data'].shape:
                        p2 = hm['p2'][g, p]
                    else:
                        p2 = hm['p2']

                    if hm['p3'].shape == hm['data'].shape:
                        p3 = hm['p3'][g, p]
                    else:
                        p3 = hm['p3']

                    # Switch to always recalculating scaled factor
                    # for each layer.
                    # Necessary for twomoment anyway
                    N = radsim.dsd(hm['data'][g, p], D, hm['type'],
                                   rho_a[g, p],
                                   self.met['Temperature'][g, p]-273.15,
                                   hm['dmin'], hm['dmax'],
                                   p1, p2, p3, hm['fc'],
                                   0, hm['apm'], hm['bpm'])
                    ze, zr, kr = radsim.zeff(self.settings['freq'],
                                             D,
                                             N,
                                             self.settings['k2'],
                                             self.settings['do_ray'],
                                             qext,
                                             qbsca)
                    # Add up the reflectivity factors
                    z_vol[g, p] += ze
                    z_ray[g, p] += zr
                    kr_vol[g, p] += kr

        Ze_ray = 10*np.log10(z_ray)
        Ze_ray[z_vol <= 0] = -999.
        Ze_non = 10*np.log10(z_vol)
        Ze_non[z_vol <= 0] = -999.

        # Calculate attenuation (gases and hydrometeors)
        a_to_vol = np.zeros(kr_vol.shape)
        g_to_vol = np.zeros(kr_vol.shape)
        g_vol = np.zeros(kr_vol.shape)
        for g in range(z_vol.shape[0]):
            for p in range(z_vol.shape[1]):
                g_vol[g, p] = radsim.gases(self.met['Pressure'][g, p],
                                           self.met['Temperature'][g, p],
                                           self.met['RH'][g, p],
                                           self.settings['freq'])
                a_to_vol[g, p] = 2 * radsim.math_lib.path_integral(
                    kr_vol[:, p],
                    self.met['Height'][:, p], 1, g)
                g_to_vol[g, p] = radsim.math_lib.path_integral(
                    g_vol[:, p],
                    self.met['Height'][:, p], 1, g)

        # Return Z effective, Z rayleigh and Z corrected.
        dBZe = Ze_non - a_to_vol - g_to_vol
        dBZe[z_vol <= 0] = -999.
        return {'Zeff': Ze_non,
                'Zray': Ze_ray,
                'Zcorr': dBZe}

    def read_mie_table(self):
        # Read Mie table
        with open(self.settings['mie_table_name'], 'rb') as file:
            data = file.read()
            mt = {}

            mt_nfreq, mt_ntt, mt_nf, mt_nd = struct.unpack('i'*4, data[4:20])

            start = 28
            inc = 8*mt_nfreq
            mt['freq'] = np.array(struct.unpack('d'*mt_nfreq,
                                                data[start:start+inc]))

            start = start + inc
            inc = 8*mt_ntt
            mt['tt'] = np.array(struct.unpack('d'*mt_ntt,
                                              data[start:start+inc]))

            start = start + inc
            inc = 8*mt_nf
            mt['f'] = np.array(struct.unpack('d'*mt_nf, data[start:start+inc]))
            mt['f_bins'] = (mt['f'][1:]+mt['f'][:-1])/2.

            start = start + inc
            inc = 4*mt_ntt
            mt['phase'] = np.array(struct.unpack('i'*mt_ntt,
                                                 data[start:start+inc]))

            start = start + inc
            inc = 8*mt_nd
            mt['D'] = np.array(struct.unpack('d'*mt_nd,
                                             data[start:start+inc]))

            start = start + inc
            inc = 8*mt_nd*mt_ntt*mt_nf*mt_nfreq
            qext = struct.unpack('d'*mt_nd*mt_ntt*mt_nf*mt_nfreq,
                                 data[start:start+inc])
            mt['qext'] = np.array(qext).reshape((mt_nfreq,
                                                 mt_nf,
                                                 mt_ntt,
                                                 mt_nd))

            start = start + inc
            inc = 8*mt_nd*mt_ntt*mt_nf*mt_nfreq
            qbsca = struct.unpack('d'*mt_nd*mt_ntt*mt_nf*mt_nfreq,
                                  data[start:start+inc])
            mt['qbsca'] = np.array(qbsca).reshape((mt_nfreq,
                                                   mt_nf,
                                                   mt_ntt,
                                                   mt_nd))

            mt_ttl = mt['tt'][mt['phase'] == 0]
            mt['liq_bins'] = (mt_ttl[1:] + mt_ttl[:-1])/2.
            mt_tti = mt['tt'][mt['phase'] == 1]
            mt['ice_bins'] = (mt_tti[1:] + mt_tti[:-1])/2.

            mt['cnt_liq'] = (mt['phase'] == 0).sum()
            mt['cnt_ice'] = (mt['phase'] == 1).sum()
            self.mt = mt

    def read_hclass(self):
        # This sets the size distributions for the different hydrometeors
        self.hclass = []
        classinput = np.genfromtxt(self.settings['hclass_file'],
                                   dtype=['f']*(len(hclasscols)-1)+['S30'])
        for i in range(len(classinput)):
            tempdict = {hclasscols[j]: classinput[i][j]
                        for j in range(len(classinput[0]))}
            if (tempdict['rho']>0) and (tempdict['apm']<0):
                tempdict['apm'] = (np.pi/6.)*tempdict['rho']
                tempdict['bpm'] = 3.
            del(tempdict['rho'])
            self.hclass.append(tempdict)

    def read_hydrodata(self):
        # Some sample hydrometeor/met data
        hm_input = np.genfromtxt(self.settings['hydro_file'], skip_header=2)
        self.met = {}
        # Setup the met data
        for i in range(4):
            self.met[metnames[i]] = hm_input[:, i][:, None]
        # Add the hydrometeor data to the hclass data
        for i in range(hm_input.shape[1]-4):
            self.hclass[i]['data'] = hm_input[:, i+4][:, None]
        self.met['Height'] /= 1000
