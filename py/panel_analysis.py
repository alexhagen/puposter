import sys
import numpy as np
from colour import Color
from dateutil import parser
import matplotlib.dates as md
from os.path import expanduser
sys.path.append(expanduser("~") + "/code")
from pymf import ctmfd as ahi
from pym import func as ahf
from pyg import twod as ahp


panel_dd = ahi.ctmfd_data()
panel_dd.add_data(expanduser("~") + "/data/tena/m41_panel_dd_1E8_200_cm_02_18_16.lvm")
panel_dd.add_data(expanduser("~") + "/data/tena/m16_panel_dd_1E8_200_cm_02_18_16.lvm")
panel_dd.add_data(expanduser("~") + "/data/tena/m42_panel_dd_1E8_200_cm_02_18_16.lvm")
panel_dd.remove_pressure(4.00)
panel_dd_curve = ahf.curve(panel_dd.p, panel_dd.wt, name="Panel - DD",
                           u_x=panel_dd.p_sigma, u_y=panel_dd.wt_sigma)

panel_dd_cf = ahi.ctmfd_data()
panel_dd_cf.add_data(expanduser("~") + "/data/tena/m41_panel_dd_1E8_200_cm_cf_100_cm_02_18_16.lvm")
panel_dd_cf.add_data(expanduser("~") + "/data/tena/m16_panel_dd_1E8_200_cm_cf_100_cm_02_18_16.lvm")
panel_dd_cf.add_data(expanduser("~") + "/data/tena/m42_panel_dd_1E8_200_cm_cf_100_cm_02_18_16.lvm")
panel_dd_cf.remove_pressure(4.00)
panel_dd_cf_curve = ahf.curve(panel_dd_cf.p, panel_dd_cf.wt,
                              name="Panel - DD + Cf", u_x=panel_dd_cf.p_sigma,
                              u_y=panel_dd_cf.wt_sigma)

I_dd = 1.0E8
I_cf = 6.0E4
rad_m16 = 1.2
l_m16 = 3.4
d_m16_cf = 69.0
d_m16_dd = 100.0
omega_m16_cf = 6.0 * rad_m16 * l_m16 / (4.0 * 3.14 * d_m16_cf**2.0)
omega_m16_dd = 6.0 * rad_m16 * l_m16 / (4.0 * 3.14 * d_m16_dd**2.0)

panel_dd_cf_eff_curve = (1.0 / panel_dd_cf_curve) / (I_cf * omega_m16_cf)
panel_dd_eff_curve = (1.0 / panel_dd_curve) / (I_dd * omega_m16_dd)

plot = panel_dd_eff_curve.plot(linecolor="#746C66", linestyle='-')
plot = panel_dd_cf_eff_curve.plot(linecolor="#E3AE24", linestyle='-', addto=plot)

plot.lines_on()
plot.xlabel('Negative Pressure ($p_{neg}$) [$bar$]')
plot.ylabel('Efficiency ($\eta$) [ ]')
plot.legend()
plot.ax.set_yscale('log')

plot.export("../img/panel_analysis", sizes=['cs'], formats=['pdf', 'pgf'],
            customsize=[4.5, 3])
