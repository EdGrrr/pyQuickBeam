import numpy as np
import struct
import radsim

metnames = ['Height', 'Pressure',
            'Temperature', 'RH']
hydronames = ['Cloud Water', 'Rain',
              'Snow', 'Aggregate',
              'Graupel', 'Ice']
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
            ('hydro_file', 'data/hydromet_example.dat'),
            ('hclass_file', 'data/hclass.dat'),
            ('mie_table_name', 'data/mie_table.dat')]

# Constants
maxhclass = 50  # Number of hm classes
nd = 85  # Number of diameters
dmin, dmax = 0.1, 10000.  # Min and max particle size
D = np.exp(np.linspace(np.log(dmin), np.log(dmax), 85))


class Quickbeam:
    def __init__(self):
        '''Populate the arrays from example files'''
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
            # Rho_eff calcualtion should only be done once
            if (hm['phase'] == 1) & (hm['rho'] == -1):
                rho_eff = (6/np.pi)*hm['apm']*(D*1e-6)**(hm['bpm']-3)
                rho_eff[rho_eff < 5] = 5.
                rho_eff[rho_eff > 917] = 917.
                hm['f'] = np.digitize(rho_eff/917., self.mt['f_bins'])
            else:
                rho_eff = np.zeros(D.shape) + hm['rho']

            for g in range(z_vol.shape[0]):
                for p in range(z_vol.shape[1]):
                    if hm['data'][g, p] == 0:
                        continue
                    if self.settings['use_mie_table'] == 1:
                        if hm['phase'] == 0:
                            temp_ind = int((self.met['Temperature'][g, p]-210.65)/5.)
                            if temp_ind < 0: temp_ind = 0
                            if temp_ind > (self.mt['cnt_liq']-1): 
                                temp_ind = self.mt['cnt_liq']-1
                            qext = self.mt['qext'][freq_ind, 0, temp_ind]
                            qbsca = self.mt['qbsca'][freq_ind, 0, temp_ind]
                        else:
                            temp_ind = int((self.met['Temperature'][g, p]
                                            -180.65)/5.)
                            if temp_ind <0: temp_ind = 0
                            if temp_ind > (self.mt['cnt_ice']-1): 
                                temp_ind = self.mt['cnt_ice']-1
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
                        qext, qbsca = radsim.optical_sphere(
                            self.settings['freq'],
                            D,
                            self.met['Temperature'][g, p]-273.15,
                            hm['phase'],
                            rho_eff)

                    if 'fc' not in hm.keys():
                        hm['fc'] = np.zeros(D.shape)
                        hm['scaled'] = 0

                    N = radsim.dsd(hm['data'][g, p], D, hm['type'],
                                   rho_a[g, p],
                                   self.met['Temperature'][g, p]-273.15,
                                   hm['dmin'],
                                   hm['dmax'], hm['rho'], hm['p1'], hm['p2'],
                                   hm['p3'], hm['fc'],
                                   hm['scaled'], hm['apm'], hm['bpm'])
                    # print N
                    ze, zr, kr = radsim.zeff(self.settings['freq'],
                                             D,
                                             N,
                                             self.settings['k2'],
                                             self.settings['do_ray'],
                                             qext,
                                             qbsca)
                    # print ze, zr, kr
                    z_vol[g, p] += ze
                    z_ray[g, p] += zr
                    kr_vol[g, p] += kr

        Ze_ray = 10*np.log10(z_ray)
        Ze_ray[z_vol <= 0] = -999.
        Ze_non = 10*np.log10(z_vol)
        Ze_non[z_vol <= 0] = -999.

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

        dBZe = Ze_non - a_to_vol - g_to_vol
        dBZe[z_vol <= 0] = -999.
        return {'Zeff': Ze_non,
                'Zray': Ze_ray,
                'Zcorr': dBZe}

    def read_mie_table(self):
        # Read Mie table
        with open(self.settings['mie_table_name']) as file:
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
        self.hclass = []
        classinput = np.genfromtxt(self.settings['hclass_file'])
        for i in range(6):
            self.hclass.append({hclasscols[j]: classinput[i, j]
                                for j in range(classinput.shape[1])})
            self.hclass[i]['name'] = hydronames[i]

    def read_hydrodata(self):
        hm_input = np.genfromtxt(self.settings['hydro_file'], skip_header=2)
        self.met = {}
        # Setup the met data
        for i in range(4):
            self.met[metnames[i]] = hm_input[:, i][:, None]
        # Add the hydrometeor data to the hclass data
        for i in range(hm_input.shape[1]-4):
            self.hclass[i]['data'] = hm_input[:, i+4][:, None]
        self.met['Height'] /= 1000