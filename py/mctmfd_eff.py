import sys
import numpy as np
from colour import Color
from dateutil import parser
import matplotlib.dates as md
from os.path import expanduser
sys.path.append(expanduser("~") + "/code")
from pymf import ctmfd as ahi
from pym import func as ahf

m16_cf = ahi.ctmfd_data()
m16_cf.add_data(expanduser("~") + '/data/tena/m16_cf_69_cm_12_09_15.lvm')
m16_cf_curve = ahf.curve(m16_cf.p, m16_cf.wt, name='Cf - 12/09/15',
                          u_x=m16_cf.p_sigma, u_y=m16_cf.wt_sigma)

m16_dd = ahi.ctmfd_data()
m16_dd.add_data(expanduser("~") + '/data/tena/m16_dd_1E7_100_cm_12_08_15.lvm')
m16_dd_curve = ahf.curve(m16_dd.p, m16_dd.wt, name='DD - 12/08/15',
                         u_x=m16_dd.p_sigma, u_y=m16_dd.wt_sigma)

m16_dd2 = ahi.ctmfd_data()
m16_dd2.add_data(expanduser("~") + '/data/tena/m16_dd_1E7_100_cm_12_21_15.lvm')
m16_dd2.add_data(expanduser("~") + '/data/tena/m16_dd_1E7_100_cm_12_22_15.lvm')
m16_dd2_curve = ahf.curve(m16_dd2.p, m16_dd2.wt, name='DD - 12/22/15',
                          u_x=m16_dd2.p_sigma, u_y=m16_dd2.wt_sigma)


m16_dd3 = ahi.ctmfd_data()
m16_dd3.add_data(expanduser("~") + '/data/tena/m16_dd_1E8_100cm_01_07_16.lvm')
m16_dd3_curve = ahf.curve(m16_dd3.p, m16_dd3.wt, name='DD - 1/08/16',
                          u_x=m16_dd3.p_sigma, u_y=m16_dd3.wt_sigma)

I_dd = 1.0E7
I_cf = 6.0E4
rad_m16 = 1.2
l_m16 = 3.4
d_m16_cf = 69.0
d_m16_dd = 100.0
omega_m16_cf = 2.0 * rad_m16 * l_m16 / (4.0 * 3.14 * d_m16_cf**2.0)
omega_m16_dd = 2.0 * rad_m16 * l_m16 / (4.0 * 3.14 * d_m16_dd**2.0)

m16_cf_eff_curve = (1.0 / m16_cf_curve) / (I_cf * omega_m16_cf)
m16_dd_eff_curve = (1.0 / m16_dd2_curve) / (I_dd * omega_m16_dd)

I_cf = 7.5E4
I_dd = 1.0E8
rad_r5 = 1.8
l_r5 = 3.6
d_r5_cf = 69.0
d_r5_dd = 100.0
omega_r5_cf = 2.0 * rad_r5 * l_r5 / (4.0 * 3.14159 * d_r5_cf**2.0)
omega_r5_dd = 2.0 * rad_r5 * l_r5 / (4.0 * 3.14159 * d_r5_dd**2.0)

r5_cf = ahi.ctmfd_data()
r5_cf.add_data(expanduser("~") + '/data/tena/r5_cf_69_cm_09_30_15.csv')
r5_cf_curve = ahf.curve(r5_cf.p, r5_cf.wt, name='Cf',
                        u_x=r5_cf.p_sigma, u_y=r5_cf.wt_sigma)
r5_cf_eff_curve = (1.0 / r5_cf_curve) / (I_cf * omega_r5_cf)
plot = r5_cf_eff_curve.plot(linecolor='#E3AE24', linestyle='-')
plot = m16_cf_eff_curve.plot(linecolor="#E3AE24", linestyle='--', addto=plot)


r5_dd = ahi.ctmfd_data()
r5_dd.add_data(expanduser("~") + '/data/tena/r5_dd_1E7_100_cm_10_06_15_1.csv')
r5_dd.add_data(expanduser("~") + '/data/tena/r5_dd_1E7_100_cm_10_06_15_2.csv')
r5_dd.add_data(expanduser("~") + '/data/tena/r5_dd_1E7_100_cm_10_01_15.csv')
r5_dd.add_data(expanduser("~") + '/data/tena/r5_dd_1E7_100_cm_10_01_15_2.csv')
r5_dd_curve = ahf.curve(r5_dd.p, r5_dd.wt, name='DD',
                        u_x=r5_dd.p_sigma, u_y=r5_dd.wt_sigma)
r5_dd_eff_curve = (1.0 / r5_dd_curve) / (I_dd * omega_r5_cf)
plot = r5_dd_eff_curve.plot(addto=plot, linecolor='#746C66', linestyle='-')
plot = m16_dd_eff_curve.plot(addto=plot, linecolor='#746C66', linestyle='--')

plot.add_data_pointer(2.75, curve=r5_cf_eff_curve, string='Cf',
                      place=(3.2, 2.0E-4))
plot.add_data_pointer(3.65, curve=r5_dd_eff_curve, string='DD',
                      place=(2.8, 3.0E-5))
plot.lines_on()
plot.ylabel('Efficiency ($\eta$) [ ]')
plot.xlabel('Pressure ($p_{neg}$) [$\mathrm{bar}$]')
plot.ax.set_yscale('log')

plot.export("../img/ctmfd_eff", sizes=['cs'], formats=['pdf', 'pgf'],
            customsize=[4.5, 3.0])
#plot.show()
plot.close()
