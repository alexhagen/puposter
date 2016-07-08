import sys
import numpy as np
from colour import Color

from os.path import expanduser
sys.path.append(expanduser("~") + "/code")
from pymf import ctmfd as ahi
from pym import func as ahm
from pyg import twod as ahp

def exp_comp(x,a,b):
    return a * (1.0 - np.exp(b * x))

v = [4.0, 15.0, 40.0, 65.0]
eff = [15.0, 36.0, 60.0, 80.0]
f = ahm.curve(v, eff, name='$\eta$')
f.fit_gen(exp_comp, guess=[100.0, -0.05])
figure = f.plot(linecolor="#746C66")
figure.lines_off()
figure = f.plot_fit(addto=figure, linecolor="#E3AE24", xmin=0.0, xmax=100.0)
figure.lines['$\eta$0'].set_linewidth(0.0)
figure.lines['$\eta$0'].set_markersize(6)
figure.lines['$\eta$fit0'].set_linewidth(1.0)
figure.lines['$\eta$fit0'].set_markersize(0)
figure.add_data_pointer(50.0, point=f.fit_at(50.0),
                        string='Interaction efficiency scales\n' +
                               'by  attenuation law',
                        place=(40.0, 20.0))
figure.ylim(0.0, 100.0)
figure.xlim(0.0, 100.0)
figure.xlabel(r'Sensitive Volume ($v$) [$cm^{3}$]')
figure.ylabel(r'Interaction Efficiency ($\varepsilon$) [$\%$]')
figure.export('../img/eff_vs_bulb_size', sizes=['cs'],
              formats=['pdf', 'pgf'], customsize=[9.0, 3.0])
