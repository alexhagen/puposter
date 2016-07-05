import sys
import numpy as np
from os.path import expanduser
sys.path.append(expanduser("~") + "/code")
from pymf import ctmfd as ahi
from pym import func as ahm
from pyg import twod as ahp

t = np.linspace(0,5.0*2.0*3.14/40000.0,num=1000)
p = np.sin(40000*t);
v = p.copy();
v[p>-0.75]=0.0;
v[p<=-0.75]=-1.0 * p[p<=-0.75];

p_curve = ahm.curve(t,p,name='$p$');
v_curve = ahm.curve(t,v,name='$V$');


figure = p_curve.plot();
figure.ax.set_yticks([-1,0,1],['$\max\,p_{neg}$','$0.0$','$\min\,p_{neg}$'])
figure = v_curve.plot(addto=figure, yy=True)
figure.ax.set_yticks([0,1],['$0.0$','$\max\,V$'])
figure.fill_between(t,np.zeros(np.size(t)),
    v,fc='#E3AE24');
figure.lines_on();
figure.markers_off();
#figure.xlabel('Time ($t$)');
figure.ylabel('Pressure Amplitude '
    +'($\left| p \\right|$)');
figure.ax2.set_ylabel('Sensitive Volume '
    +'($V$)');
figure.ax.set_xticks([],[]);
figure.ax2.set_ylim(-1,1);

figure.add_subplot(212);
phi = np.zeros(t.size);
t_int = t<1.0*3.14/40000.0;
t_prompt = np.all([[t>1.0*3.14/40000.0],[t<1.1*3.14/40000.0]],axis=0).reshape((1000,));
t_fiss = np.all([[t>1.1*3.14/40000.0],[t<6.0*3.14/40000.0]],axis=0).reshape((1000,));
t_int2 =np.all([[t>8.0*3.14/40000.0],[t<9.0*3.14/40000.0]],axis=0).reshape((1000,));
t_fiss2 = np.all([[t>9.1*3.14/40000.0],[t<14.0*3.14/40000.0]],axis=0).reshape((1000,));
t_prompt2 = np.all([[t>9.0*3.14/40000.0],[t<9.1*3.14/40000.0]],axis=0).reshape((1000,));
phi[t_int] = 1.0;
phi[t_fiss] = 0.5*np.exp(-10000.0*(t[t_fiss]-1.0*3.14/40000.0));
phi[t_int2] = 1.0;
phi[t_prompt] = 2.0;
phi[t_prompt2] = 2.0;
phi[t_fiss2] = 0.5*np.exp(-10000.0*(t[t_fiss2]-9.0*3.14/40000.0));
mid = np.zeros(t.size);
mid[t_int] = 0.5;
mid[t_fiss] = 0.25*np.exp(-10000.0*(t[t_fiss]-1.0*3.14/40000.0));
mid[t_int2] = 0.5;
mid[t_prompt] = 1.0;
mid[t_prompt2] = 1.0;
mid[t_fiss2] = 0.25*np.exp(-10000.0*(t[t_fiss2]-9.0*3.14/40000.0))
phi_curve = ahm.curve(t,phi,name='$\phi$')
mid_curve = ahm.curve(t,mid,name='mid')

figure = phi_curve.plot(addto=figure, axes=figure.ax_subp[0]);
figure.fill_between(t[t_int],
    np.zeros(np.size(t[t_int])),
    phi[t_int],fc='#A7A9AC',axes=figure.ax_subp[0]);
figure.fill_between(t[t_prompt],
    np.zeros(np.size(t[t_prompt])),
    phi[t_prompt],fc='#2EAFA4',axes=figure.ax_subp[0]);
figure.fill_between(t[t_prompt2],
    np.zeros(np.size(t[t_prompt2])),
    phi[t_prompt2],fc='#2EAFA4',axes=figure.ax_subp[0]);
figure.fill_between(t[t_fiss],
    np.zeros(np.size(t[t_fiss])),
    phi[t_fiss],fc='#5C8727',axes=figure.ax_subp[0]);
figure.fill_between(t[t_int2],
    np.zeros(np.size(t[t_int2])),
    phi[t_int2],fc='#A7A9AC',axes=figure.ax_subp[0]);
figure.fill_between(t[t_fiss2],
    np.zeros(np.size(t[t_fiss2])),
    phi[t_fiss2],fc='#5C8727',axes=figure.ax_subp[0]);
figure.add_data_pointer(0.5*3.14/40000.0,mid_curve,
    'Interrogation Pulse',place=(1.5*3.14/40000.0,0.75),
    axes=figure.ax_subp[0]);
figure.add_data_pointer(1.05*3.14/40000.0,mid_curve,
    'Prompt Signal',place=(2.0*3.14/40000.0,1.0),
    axes=figure.ax_subp[0]);
figure.add_data_pointer(2.5*3.14/40000.0,mid_curve,
    'Die Away Signal',place=(4.5*3.14/40000.0,0.5),
    axes=figure.ax_subp[0])
figure.xlabel('Time ($t$)', figure.ax_subp[0])
figure.ylabel('Flux ($\phi$)', figure.ax_subp[0])
figure.ax_subp[0].set_xticks([], [])
figure.ax_subp[0].set_yticks([0.25, 1.0], ['$\phi_{f}$', '$\phi_{int}$'])
figure.ax_subp[0].set_ylim(0, 1.25)

figure.markers_off()
figure.lines_on()
figure.export('../img/atmfd_p', sizes=['cs'],
              formats=['pdf', 'pgf'], customsize=[9.0, 6.0])
