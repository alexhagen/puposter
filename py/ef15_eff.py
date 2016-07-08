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


pwr_new, eff_new, u_eff = \
    np.loadtxt('Fig7-Cf-Eff.csv',
               delimiter=",", dtype=object,
               converters={0: float, 1: float, 2: float},
               unpack=True)
# points = [range(0, 6), range(7, 11), range(13, 17), range(19, 25),
#           range(26, 31), range(32, 39), range(40, 48), range(49, 57),
#           range(58, 63), range(64, 68), range(69, 82)]
# pwr = [np.mean(pwr_new[idx]) for idx in points]
# u_pwr = [np.std(pwr_new[idx]) for idx in points]
# eff = np.array([np.mean(eff_new[idx]) for idx in points])
# u_eff = np.array([np.mean(u_eff[idx]) for idx in points])
# std_eff = np.array([np.std(eff_new[idx]) for idx in points])
# u_eff = np.sqrt(np.power(u_eff, 2.0) + np.power(std_eff, 2.0))
pwr = pwr_new
eff = eff_new
cf_tena_eta = ahf.curve(pwr, 100. * eff, u_y=u_eff, name="cf_tena_setup")

pwr_new, eff_new, u_eff = \
    np.loadtxt('Fig7-DD-Eff.csv',
               delimiter=",", dtype=object,
               converters={0: float, 1: float, 2: float},
               unpack=True)
# pwr = [np.mean(pwr_new[0:1]), np.mean(pwr_new[2]), np.mean(pwr_new[3]),
#        np.mean(pwr_new[4:])]
# u_pwr = [np.std(pwr_new[0:1]), np.std(pwr_new[2]), np.std(pwr_new[3]),
#          np.std(pwr_new[4:])]
# eff = np.array([np.mean(eff_new[0:1]), np.mean(eff_new[2]), np.mean(eff_new[3]),
#        np.mean(eff_new[4:])]) / nps_dd
# u_eff = np.array([np.mean(u_eff[0:1]), np.mean(u_eff[2]), np.mean(u_eff[3]),
#          np.mean(u_eff[4:])]) / nps_dd
# std_eff = np.array([np.std(eff_new[0:1]), np.std(eff_new[2]), np.std(eff_new[3]),
#            np.std(eff_new[4:])]) / nps_dd
# u_eff = np.sqrt(np.power(u_eff, 2.0) + np.power(std_eff, 2.0))
pwr = pwr_new
eff = eff_new
dd_tena_eta = ahf.curve(pwr, 100. * eff, u_y=u_eff, name="dd_tena_setup")

plot = ahp.ah2d()
plot = dd_tena_eta.plot(linecolor='#746C66', linestyle='-')
plot = cf_tena_eta.plot(linecolor='#E3AE24', linestyle='-', addto=plot)
plot.lines_off()
plot.ax.set_yscale('log')
plot.add_data_pointer(3.5, curve=cf_tena_eta, string=r"$Cf$",
                      place=(2.2, 2E-2))
plot.add_data_pointer(6.5, curve=dd_tena_eta, string=r"$DD$",
                      place=(7.5, 1E-5))
plot.xlabel(r'Power ($P$) [$W$]')
plot.ylabel(r'Efficiency ($\eta$) [$\%$]')
plot.export('../img/ef15_eff', formats=['pdf', 'pgf'], sizes=['cs'],
            customsize=(4.5, 3.0))
