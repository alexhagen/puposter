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
