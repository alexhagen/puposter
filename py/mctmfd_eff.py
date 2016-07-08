import sys
import numpy as np
from colour import Color
from dateutil import parser
import matplotlib.dates as md
from os.path import expanduser
sys.path.append(expanduser("~") + "/code")
from pymf import ctmfd as ahi
from pym import func as ahf

pneg, eff, u_eff = \
    np.loadtxt('Fig8-Cf-Eff.csv',
               delimiter=",", dtype=object,
               converters={0: float, 1: float, 2: float},
               unpack=True)
m16_cf_curve = ahf.curve(pneg, 100. * eff, u_y=u_eff, name='$Cf$')
pneg, eff, u_eff = \
    np.loadtxt('Fig8-DD-Eff.csv',
               delimiter=",", dtype=object,
               converters={0: float, 1: float, 2: float},
               unpack=True)
m16_dd_curve = ahf.curve(pneg, 100. * eff, u_y=u_eff, name='$DD$')
pneg, r, u_r = \
    np.loadtxt('Fig8-Ratio.csv',
               delimiter=",", dtype=object,
               converters={0: float, 1: float, 2: float},
               unpack=True)
m16_ratio_curve = ahf.curve(pneg, r, u_y=u_r, name='$r$')

plot = m16_cf_curve.plot(linestyle='-', linecolor='#E3AE24')
plot = m16_dd_curve.plot(linestyle='-', linecolor='#746C66', addto=plot)
plot.add_data_pointer(3.0, curve=m16_cf_curve, string='Cf',
                      place=(3.2, 2.0E-2))
plot.add_data_pointer(3.2, curve=m16_dd_curve, string='DD',
                      place=(2.8, 3.0E-3))
plot.lines_on()
plot.ylabel('Efficiency ($\eta$) [ ]')
plot.xlabel('Pressure ($p_{neg}$) [$\mathrm{bar}$]')
plot.ax.set_yscale('log')
plot.xlim(2.6, 4.2)
plot.ylim(1.E-5, 5.E0)

plot.export("../img/m16_eff", sizes=['cs'], formats=['pdf', 'pgf'],
            customsize=[4.5, 3.0])

pneg, u_pneg, eff, u_eff = \
    np.loadtxt('Fig9-Cf-Eff.csv',
               delimiter=",", dtype=object,
               converters={0: float, 1: float, 2: float, 3: float},
               unpack=True)
r5_cf_curve = ahf.curve(pneg, 100. * eff, u_x=u_pneg, u_y=u_eff, name='$Cf$')
pneg, u_pneg, eff, u_eff = \
    np.loadtxt('Fig9-DD-Eff.csv',
               delimiter=",", dtype=object,
               converters={0: float, 1: float, 2: float, 3: float},
               unpack=True)
r5_dd_curve = ahf.curve(pneg, 100. * eff, u_x=u_pneg, u_y=u_eff, name='$DD$')
pneg, u_pneg, r, u_r = \
    np.loadtxt('Fig9-Ratio.csv',
               delimiter=",", dtype=object,
               converters={0: float, 1: float, 2: float, 3: float},
               unpack=True)
r5_ratio_curve = ahf.curve(pneg, r, u_x=u_pneg, u_y=u_r, name='$r$')

plot = r5_cf_curve.plot(linestyle='-', linecolor='#E3AE24')
plot = r5_dd_curve.plot(linestyle='-', linecolor='#746C66', addto=plot)
plot.add_data_pointer(2.75, curve=r5_cf_curve, string='Cf',
                      place=(3.2, 2.0E-2))
plot.add_data_pointer(3.65, curve=r5_dd_curve, string='DD',
                      place=(2.8, 3.0E-3))
plot.lines_on()
plot.ylabel('Efficiency ($\eta$) [$\%$]')
plot.xlabel('Pressure ($p_{neg}$) [$\mathrm{bar}$]')
plot.ax.set_yscale('log')
plot.xlim(2.6, 4.2)
plot.ylim(1.E-5, 5.E0)

plot.export("../img/r5_eff", sizes=['cs'], formats=['pdf', 'pgf'],
            customsize=[4.5, 3.0])
