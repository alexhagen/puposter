import sys
import numpy as np
from os.path import expanduser
sys.path.append(expanduser("~") + "/code")
from pymf import ctmfd as ahsoa
from pym import func as ahm

E_high = 6.00;

p = [5.0,6.0,7.0,10.0];
wt = [186.31,4.53,2.81,0.00];
wt_err = [0.00,4.53,1.99,0.00];
p_u = [5.0];
wt_u = [46.216];
wt_u_err = [20.668];
p_be = [5.0];
wt_be = [6.882];
wt_be_err = [3.07];
detection = ahm.curve(p,wt,name='connector',u_y=wt_err);
detection_plot = detection.plot(linecolor='#746C66', linestyle='-', legend=False);
u_det = ahm.curve(p_u,wt_u,name='U-235',u_y=wt_u_err);
be_det = ahm.curve(p_be,wt_be,name='Be',u_y=wt_be_err);
detection_plot = u_det.plot(linecolor='#E3AE24', addto=detection_plot);
detection_plot = be_det.plot(linecolor='#2EAF9B', addto=detection_plot);
detection_plot.add_text(7.25, 100.0, string="Detection of \n photoneutron contamination")
detection_plot.markers_on();
detection_plot.lines_off();
detection_plot.lines['connector'].set_linewidth(1.0);
detection_plot.lines['connector'].set_markersize(0)
detection_plot.ylim(0,200.0);
detection_plot.xlim(0,10.0);
detection_plot.add_data_pointer(p_u[0],point=wt_u[0],
	string="NU Detection",place=(3.0,100.0));
detection_plot.add_data_pointer(p_be[0],point=wt_be[0],
	string="Be Detection",place=(3.5,30.0));
detection_plot.fill_between(p, wt, 200.0 * np.ones_like(wt), fc="#D1D3D4", ec="#D1D3D4", name="control")
detection_plot.xlabel('Pressure ($p$) [$bar$]');
detection_plot.ylabel('Waiting Time ($t_{wait}$) [$s$]');
detection_plot.legend();
detection_plot.export('../img/photoneutron_detection',sizes=['cs'],
	formats=['pdf', 'pgf'],customsize=[9.0, 3.0]);
