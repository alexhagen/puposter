import sys
import numpy as np
from colour import Color
from dateutil import parser
from datetime import datetime
from matplotlib.dates import strpdate2num
from os.path import expanduser
sys.path.append(expanduser("~") + "/code")
from pymf import ctmfd as ahi
from pym import func as ahf
from pyg import twod as ahp
from scipy.cluster.vq import kmeans, vq, whiten
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import re

class datmfd_analysis():
    def __init__(self, fnames, sources, name, nps, u_nps):
        self.export_name = name
        # first check all of the filenames and cluster the power and detections
        levels = self.cluster_power(fnames, nps, u_nps)
        self.p = [x.p for x in levels]
        self.u_p = [x.u_p for x in levels]
        self.cr = [x.cr for x in levels]
        self.u_cr = [x.u_cr for x in levels]
        # then, we want to calculate the efficiency
        for crs in self.cr:
            self.eta = np.array(self.cr) / nps
            self.u_eta = np.array(self.eta) * \
                np.sqrt(np.divide(self.u_cr, self.cr)**2.0 +
                        np.divide(u_nps, nps)**2.0)

    def return_cr(self):
        return ahf.curve(self.p, self.cr, u_x=self.u_p, u_y=self.u_cr,
                         name=self.export_name + '_cr')
    def return_eff(self):
        return ahf.curve(self.p, self.eta, u_x=self.u_p, u_y=self.u_eta,
                         name=self.export_name + '_eta')

    def cluster_power(self, fnames, nps, u_nps):
        t_pwr = []
        pwr = []
        u_pwr = []
        t_det = []
        det = []
        levels = []
        # plot the average power
        plot = ahp.ah2d()
        colors = ['#E3AE24', '#B63F97', '#5C6F7B', '#5C8727', '#7E543A', '#A3792C', '#7299C6', '#B8B308']
        for fname in fnames:
            # check if an alex fileset or a brian file
            if "pwr" in fname:
                # open power file
                t_pwr_new, pwr_new, u_pwr_new = \
                    np.loadtxt(expanduser("~") + '/data/tena/' + fname + '.dat',
                               delimiter=",", dtype=object,
                               converters={0: strpdate2num("%d-%b-%Y %H:%M:%S"),
                                           1: float, 2: float}, unpack=True)
                t_pwr_new = np.array(t_pwr_new) * 3600. * 24.
                t_det_new, det_new, _ = \
                    np.loadtxt(expanduser("~") + "/data/tena/" + fname.replace("pwr", "det") + '.dat',
                               delimiter=",", dtype=object,
                               converters={0: strpdate2num("%d-%b-%Y %H:%M:%S"),
                                           1: float, 2: float}, unpack=True)
                t_det_new = np.array(t_det_new) * 3600. * 24.
                eps = 0.1
                min_samples = 10
            else:
                # get the date from the filename
                # find string matching %i_min_%2i_%2i_%2i[_%i] at end of filename
                # get the count interval from the filename
                # find string matching %i_min at end of string
                pattern = \
                    ".*_([0-9]*)_min_([0-9]*)_([0-9]*)_([0-9]*)(?:_([0-9]*))?$"
                result = re.split(pattern, fname)
                mm = result[2]
                dd = result[3]
                yy = result[4]
                t_int = float(result[1]) * 60.
                epoch = datetime.utcfromtimestamp(0)
                t_start = (datetime.strptime(
                    mm + "/" + dd + "/" + yy + " 01:00:00",
                    "%m/%d/%y %H:%M:%S") - epoch).total_seconds()
                pwr_new, eff_new, u_eff_new = \
                    np.loadtxt(expanduser("~") + '/data/tena/' + fname + '.dat',
                               delimiter=",", dtype=object,
                               converters={0: float, 1: float, 2: float},
                               unpack=True)
                u_pwr_new = 0.2 * np.ones_like(pwr_new)
                det_new = nps * np.array(eff_new)
                u_eff_new = u_nps * np.array(u_eff_new)
                t = t_start
                t_pwr_new = []
                for i in range(len(pwr_new)):
                    t_pwr_new = np.append(t_pwr_new, t)
                    t += t_int
                t_det_new = t_pwr_new.copy()
                eps = 0.24
                min_samples = 2

            # put the data in row form
            pwr_data = \
                [[float(t_pwr_new[i]), float(pwr_new[i])] for i in range(len(pwr_new))]

            # use the dbscan clustering algorithm
            X = StandardScaler().fit_transform(pwr_data)
            db = DBSCAN(eps=eps, min_samples=min_samples).fit(X)
            core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True
            labels = db.labels_

            # Number of clusters in labels, ignoring noise if present.
            n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
            marked = np.zeros_like(t_pwr_new)
            for i in range(len(t_pwr_new))[::-1]:
                label = labels[i]
                other_labels = [x for x in np.unique(labels) if x != label]
                for other in other_labels:
                    if (other >= 0) and inrange(t_pwr_new[i], t_pwr_new[labels == other]):
                        marked[i] = 1
            labels[marked == 1] = -1
            for i in [x for x in np.unique(labels) if x != -1]:
                if len(pwr_new[labels == i]) >= min_samples:
                    level = power_level(t_pwr_new[labels == i], pwr_new[labels == i],
                                        t_det_new, det_new)
                    levels = np.append(levels, level)
                    t_det = np.append(t_det, t_det_new)
                    t_pwr = np.append(t_pwr, t_pwr_new)
                    pwr = np.append(pwr, pwr_new)
                    u_pwr = np.append(u_pwr, u_pwr_new)
                    det = np.append(det, det_new)

        return levels


class power_level():
    def __init__(self, t, pwr, t_det, det):
        self.t_min = np.min(t)
        self.t_max = np.max(t)
        self.p = np.mean(pwr)
        self.u_p = np.std(pwr)
        self.t_det = t_det[np.logical_and(t_det >= self.t_min, t_det <= self.t_max)]
        self.det = det[np.logical_and(t_det >= self.t_min, t_det <= self.t_max)]
        self.cr = np.sum(self.det) / (self.t_max - self.t_min)
        self.u_cr = np.sqrt(np.sum(self.det)) / (self.t_max - self.t_min)

def inrange(x, data):
    return (x >= np.min(data)) and (x <= np.max(data))

# Tena Analysis
I_dd = 1.0E8
u_I_dd = 0.2 * 1.0E8
I_cf = 6.4E4
u_I_cf = 0.2 * 6.4E4
r_sv = 1.5
u_r_sv = 0.1
h_sv = 5.0
u_h_sv = 0.1
d_dd = 100.0
u_d_dd = 0.5
d_cf = 15.0
u_d_cf = 0.5
omega_dd = 2.0 * r_sv * h_sv / (4.0 * np.pi * d_dd**2.0)
u_omega_dd = omega_dd * np.sqrt((u_r_sv / r_sv)**2.0 + (u_h_sv / h_sv)**2.0 +
                                2.0 * (u_d_dd / d_dd)**2.0)
nps_dd = I_dd * omega_dd
u_nps_dd = nps_dd * np.sqrt((u_I_dd / I_dd)**2.0 +
                            (u_omega_dd / omega_dd)**2.0)

omega_cf = 2.0 * r_sv * h_sv / (4.0 * np.pi * d_cf**2.0)
u_omega_cf = omega_cf * np.sqrt((u_r_sv / r_sv)**2.0 + (u_h_sv / h_sv)**2.0 +
                                2.0 * (u_d_cf / d_cf)**2.0)

nps_cf = I_cf * omega_cf
u_nps_cf = nps_cf * np.sqrt((u_I_cf / I_cf)**2.0 +
                              (u_omega_cf / omega_cf)**2.0)


nps_dd_cf = nps_dd + nps_cf
u_nps_dd_cf = np.sqrt((u_nps_dd)**2.0 + (u_nps_cf)**2.0)


pwr_new, eff_new, u_eff = \
    np.loadtxt(expanduser("~") + '/data/tena/' + 'ef17_cf_15_cm_04_26_16' + '.dat',
               delimiter=",", dtype=object,
               converters={0: float, 1: float, 2: float},
               unpack=True)
pwr = [np.mean(pwr_new[0:3]), np.mean(pwr_new[4:10]), np.mean(pwr_new[11:15]),
       np.mean(pwr_new[16:])]
u_pwr = [np.std(pwr_new[0:3]), np.std(pwr_new[4:10]), np.std(pwr_new[11:15]),
         np.std(pwr_new[16:])]
eff = np.array([np.mean(eff_new[0:3]), np.mean(eff_new[4:10]), np.mean(eff_new[11:15]),
       np.mean(eff_new[16:])]) / nps_cf
u_eff = np.array([np.mean(u_eff[0:3]), np.mean(u_eff[4:10]), np.mean(u_eff[11:15]),
         np.mean(u_eff[16:])]) / nps_cf
std_eff = np.array([np.std(eff_new[0:3]), np.std(eff_new[4:10]), np.std(eff_new[11:15]),
           np.std(eff_new[16:])]) / nps_cf
u_eff = np.sqrt(np.power(u_eff, 2.0) + np.power(std_eff, 2.0))
cf_tena_eta = ahf.curve(pwr, eff, u_x=u_pwr, u_y=u_eff, name="cf_tena_setup")

pwr_new, eff_new, u_eff = \
    np.loadtxt(expanduser("~") + '/data/tena/' + 'ef17_dd_1E8_100_cm_04_27_16' + '.dat',
               delimiter=",", dtype=object,
               converters={0: float, 1: float, 2: float},
               unpack=True)
pwr = [np.mean(pwr_new[0:1]), np.mean(pwr_new[2]), np.mean(pwr_new[3]),
       np.mean(pwr_new[4:])]
u_pwr = [np.std(pwr_new[0:1]), np.std(pwr_new[2]), np.std(pwr_new[3]),
         np.std(pwr_new[4:])]
eff = np.array([np.mean(eff_new[0:1]), np.mean(eff_new[2]), np.mean(eff_new[3]),
       np.mean(eff_new[4:])]) / nps_dd
u_eff = np.array([np.mean(u_eff[0:1]), np.mean(u_eff[2]), np.mean(u_eff[3]),
         np.mean(u_eff[4:])]) / nps_dd
std_eff = np.array([np.std(eff_new[0:1]), np.std(eff_new[2]), np.std(eff_new[3]),
           np.std(eff_new[4:])]) / nps_dd
u_eff = np.sqrt(np.power(u_eff, 2.0) + np.power(std_eff, 2.0))
dd_tena_eta = ahf.curve(pwr, eff, u_x=u_pwr, u_y=u_eff, name="dd_tena_setup")

plot = dd_tena_eta.plot(linecolor='#746C66', linestyle='-')
plot = cf_tena_eta.plot(linecolor='#E3AE24', linestyle='-', addto=plot)
plot.lines_on()
plot.ax.set_yscale('log')
plot.add_data_pointer(4.0, curve=cf_tena_eta, string=r"$Cf$",
                      place=(2.7, 2E-4))
plot.add_data_pointer(6.5, curve=dd_tena_eta, string=r"$DD$",
                      place=(7.5, 1E-7))
plot.xlabel(r'Power ($P$) [$W$]')
plot.ylabel(r'Efficiency ($\eta$) [ ]')
plot.export('../img/ef15_eff', formats=['pdf', 'pgf'], sizes=['cs'],
            customsize=(6, 3.0))
plot.show()

'''
import sys
import numpy as np
from os.path import expanduser
sys.path.append(expanduser("~") + "/code")
from pymf import ctmfd as ahi
from pym import func as ahf

arr = np.loadtxt(expanduser("~") + '/data/tena/ef15_cf_69cm_12_10_15_win.csv', delimiter=',',
                 usecols=(0, 1, 2, 3, 4))
cf_pwr = arr[:, 0]
cf_cnt = arr[:, 1]
cf_u_cnt = arr[:, 2]
cf_eff = arr[:, 3]
cf_u_eff = arr[:, 4]

cf = ahf.curve(cf_pwr, cf_eff, u_y=cf_u_eff, name='Cf - 69 cm')

plot = cf.plot(linecolor="#E3AE24", linestyle="-")

arr = np.loadtxt(expanduser("~") + '/data/tena/ef15_cf_15cm_12_10_15_win.csv', delimiter=',',
                 usecols=(0, 1, 2, 3, 4))
cf_15_pwr = arr[:, 0]
cf_15_cnt = arr[:, 1]
cf_15_u_cnt = arr[:, 2]
cf_15_eff = arr[:, 3]
cf_15_u_eff = arr[:, 4]

cf_15 = ahf.curve(cf_15_pwr, cf_15_eff, u_y=cf_15_u_eff, name='Cf - 15 cm')

plot = cf_15.plot(addto=plot, linecolor="#E3AE24", linestyle="-")

arr = np.loadtxt(expanduser("~") + '/data/tena/ef15_dd_1E7_100cm_12_10_15_win.csv', delimiter=',',
                 usecols=(0, 1, 2, 3, 4))
dd_pwr = arr[:, 0]
dd_cnt = arr[:, 1]
dd_u_cnt = arr[:, 2]
dd_eff = arr[:, 3]
dd_u_eff = arr[:, 4]

dd = ahf.curve(dd_pwr, dd_eff, u_y=dd_u_eff,
               name=r"DD ($10^{7}\,\frac{n}{s}$) - 100 cm")

plot = dd.plot(addto=plot, linecolor="#746C66", linestyle="-")

plot.legend()
plot.ax.set_yscale("log", nonposy='clip')
plot.ylim(10.0**-5, 0.01)
plot.xlabel('Power Level ($P$) [$W$]')
plot.ylabel(r"Intrinsic Detection Efficiency ($\eta$) [ ]")
plot.export("../img/ef15_eff", sizes=['cs'], formats=['pdf', 'pgf'],
            customsize=[6.0, 3.0])
'''
