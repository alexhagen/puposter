import sys
import numpy as np
from os.path import expanduser
sys.path.append(expanduser("~") + "/code")
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pymf import ctmfd as ahi
from pym import func as ahm
from pyg import twod as ahp
from mpl_toolkits.axes_grid.inset_locator import inset_axes

nps = 4756271.0

e_132 = np.array([1.0000E-08, 2.0000E-01, 4.0000E-01, 6.0000E-01, 8.0000E-01,
                  1.0000E+00, 1.2000E+00, 1.4000E+00, 1.6000E+00, 1.8000E+00,
                  2.0000E+00, 2.2000E+00, 2.4000E+00, 2.6000E+00, 2.8000E+00])
n_132 = np.array([8.92392E-08, 3.21572E-06, 1.99213E-07, 1.91666E-07, 1.31752E-07,
                  1.53087E-07, 6.39654E-08, 3.19446E-08, 3.29229E-08, 4.86987E-08,
                  0.00000E+00, 6.65073E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00])
u_n_132 = np.array([0.5187, 0.0958, 0.3564, 0.3921, 0.4331, 0.4774, 0.5068, 0.8091,
                    1.0000, 0.6180, 0.0000, 1.0000, 0.0000, 0.0000, 0.0000])

e_surf = np.array([1.0000E-04, 2.0010E-01, 4.0010E-01, 6.0010E-01, 8.0010E-01,
                   1.0001E+00, 1.2001E+00, 1.4001E+00, 1.6001E+00, 1.8001E+00,
                   2.0001E+00, 2.2001E+00, 2.4001E+00, 2.6001E+00, 2.8001E+00])
n_surf = np.array([1.94447E-01, 2.60177E-01, 4.27055E-02, 2.45054E-02,
                   2.43836E-02, 1.30397E-02, 9.24296E-03, 8.60044E-03,
                   7.69931E-03, 5.57538E-03, 3.23034E-03, 1.68375E-03,
                   4.15662E-04, 2.96451E-05, 0.00000E+00]) / (100.0 * 100.0)
u_n_surf = np.array([0.0009, 0.0008, 0.0022, 0.0029, 0.0029, 0.0040, 0.0047,
                     0.0049, 0.0052, 0.0061, 0.0081, 0.0112, 0.0225, 0.0848,
                     0.0000]) / (100.0 * 100.0)

e_0 = np.array([1.0000E-08, 2.0000E-01, 4.0000E-01, 6.0000E-01, 8.0000E-01,
                1.0000E+00, 1.2000E+00, 1.4000E+00, 1.6000E+00, 1.8000E+00,
                2.0000E+00, 2.2000E+00, 2.4000E+00, 2.6000E+00, 2.8000E+00])
n_0 = np.array([1.04218E-07, 1.42173E-05, 1.36260E-06, 7.63546E-07,
                6.47598E-07, 1.71009E-07, 2.16882E-07, 1.57995E-07,
                1.16909E-07, 6.83332E-08, 6.82942E-08, 3.41746E-08,
                0.00000E+00, 0.00000E+00, 0.00000E+00])
u_n_0 = np.array([0.4900, 0.0476, 0.1548, 0.2102, 0.2265, 0.4472, 0.3899,
                  0.4534, 0.5222, 0.7071, 0.7071, 1.0000, 0.0000, 0.0000,
                  0.0000])

e_52 = np.array([1.0000E-08, 2.0000E-01, 4.0000E-01, 6.0000E-01, 8.0000E-01,
                 1.0000E+00, 1.2000E+00, 1.4000E+00, 1.6000E+00, 1.8000E+00,
                 2.0000E+00, 2.2000E+00, 2.4000E+00, 2.6000E+00, 2.8000E+00])
n_52 = np.array([4.91586E-08, 1.06091E-05, 1.11736E-06, 5.39662E-07,
                 5.56099E-07, 1.77157E-07, 1.59955E-07, 1.79698E-07,
                 1.32121E-07, 8.01516E-08, 3.85630E-08, 6.92248E-08,
                 0.00000E+00, 0.00000E+00, 0.00000E+00])
u_n_52 = np.array([0.5288, 0.0538, 0.1673, 0.2384, 0.2289, 0.4050, 0.4480,
                   0.4208, 0.5182, 0.5116, 1.0000, 0.6735, 0.0000, 0.0000,
                   0.0000])

r = [0., 52., 132.]
thresh = 0.250 # estimating the threshold to be around 500 kev
n = np.array([np.sum(n_0[np.where(e_0 > thresh)]),
              np.sum(n_52[np.where(e_52 > thresh)]),
              np.sum(n_132[np.where(e_132 > thresh)])]) / nps
print n
mcnp = ahm.curve(r, n, name='simulation')

d = [0.0, 52.0, 100.0, 132.0]
wt = np.array([2.75, 56.69, 277.15, 3.0 * 121.45])
u_wt = np.array([1.59, 28.34, 277.15, 121.45]) / wt
cr = 1.0 / wt
u_cr = cr * u_wt
u_cr_up = u_cr
u_cr_down = u_cr.copy()
u_cr_down[3] = np.inf
u_cr = np.vstack((u_cr_down, u_cr_up))

expt = ahm.curve(d, cr, u_y=u_cr, name='experiment')

x = np.arange(1.0, 140.)
y = cr[0] / np.power(x, 1.0)
r2 = ahm.curve(x, y, name=r'$\frac{1}{r}$')
plot = r2.plot(linecolor='#A7A9AC', linestyle='-')
plot = expt.plot(linecolor='#E3AE24', linestyle='-', addto=plot)
print plot.lines

#plot.add_reg_line(d, cr, yerr=u_cr, regtype='2')
# plot = mcnp.plot(linecolor='#2EAFA4', linestyle='-', yy=True, addto=plot)

plot.xlabel(r'Distance from CLINAC Axis ($d$) [$\mathrm{cm}$]')
plot.ylabel(r'Count Rate ($\dot{c}$) [$cps$]')
#plot.ylim(0.0, 0.5)
plot.xlim(-1, 140)
plot.lines_on()
plot.markers_off()
plot.lines['experiment'].set_linewidth(0.0)
plot.lines['experiment'].set_alpha(1.0)
plot.lines['experiment'].set_markersize(6)
plot.legend()
plot.ax.set_yscale("log", nonposy='clip')

n_ax = plot.ax
plot.export('n_plot', formats=['pdf'], sizes=['cs'], customsize=(4.5, 3.0))

plot = ahp.ah2d()
plot.fill_between(e_surf, np.zeros_like(e_surf), 2.E10 * 0.000115 * n_surf, fc='#A7A9AC',
                  name='floor')
plot.fill_between(e_0, np.zeros_like(e_0), 2.E10 * 0.000115 * n_0, fc="#746C66", name='on axis')
plot.fill_between(e_52, np.zeros_like(e_52), 2.E10 * 0.000115 * n_52, fc="#E3AE24",
                  name=r'$52\,cm$')
plot.fill_between(e_132, np.zeros_like(e_132), 2.E10 * 0.000115 * n_132, fc="#2EAFA4",
                  name=r'$132\,cm$')

plot.ylabel(r'Flux ($\phi$) [$\frac{n}{cm^{2}s}$]')
plot.xlabel(r'Energy ($E$) [$MeV$]')
plot.legend()

spectrum_ax = plot.ax
plot.export('spectrum_plot', formats=['pdf'], sizes=['cs'], customsize=(4.5, 3.0))


plot = ahp.ah2d()
# The concrete
plot.fill_between([-100., 100.], [-130.0, -130.0], [-100.0, -100.0],
                  fc='#D1D3D4')
plot.fill_between([-100., 100.], [-280.0, -280.0], [-130.0, -130.0],
                  fc='#C4B5AF')

with open("/Users/ahagen/code/mcnp_companion/tests/clinac_mcnp" + "/photo.out_temp", "w"):
  pass

with open("/Users/ahagen/code/mcnp_companion/tests/clinac_mcnp" + "/archive/clinac_no_target_h2_photo.out", "r") as f:
  cleanstring = f.read().replace("\n", " ").replace("_", "\n")
  with open("/Users/ahagen/code/mcnp_companion/tests/clinac_mcnp" + "/photo.out_temp", "a") as f2:
      f2.write(cleanstring)

E_p, _, zaid, _, x, y, z, _, _, _, _, _, ux, uy, uz, _, _ = \
  np.loadtxt("/Users/ahagen/code/mcnp_companion/tests/clinac_mcnp" + "/photo.out_temp", unpack=True)
xyz = np.vstack([x, y, z])
num_photoneutrons = 2.1E10 / 1.0E9 * 0.000115 * len(x)
kde = stats.gaussian_kde(xyz)

# Evaluate kde on a grid
xmin, zmin = x.min(), z.min()
xmax, zmax = x.max(), z.max()
xi = np.linspace(-100, 100, 25)
zi = np.linspace(-280, -100, 25)
X, Z = np.meshgrid(xi, zi)
ys = 2.1E10 / 1.0E9 * 0.000115 * np.array([kde((xii, 0., zii)) for xii,zii in zip(np.ravel(X), np.ravel(Z))])
Y = ys.reshape(X.shape)

cdict = {"red": ((0.0, 1.0, 1.0),
                 (1.0, 46./255., 46./255.)),
         "green": ((0.0, 1.0, 1.0),
                   (1.0, 175./255., 175./255.)),
         "blue": ((0.0, 1.0, 1.0),
                  (1.0, 164./255., 164./255.)),
         "alpha": ((0.0, 0.0, 0.0),
                   (1., 1.0, 1.0))}

fountain_run = LinearSegmentedColormap("Fountain_Run", cdict)
plt.register_cmap(cmap=fountain_run)


levels = np.linspace(np.min(ys), np.max(ys), 10)
surf = plot.ax.contourf(X, Z, Y, levels, cmap=fountain_run)

plot.add_line([-1.4, 1.4, 1.4, -1.4, -1.4], [-1.4, -1.4, 1.4, 1.4, -1.4],
              linestyle='-', linecolor='#746C66')
plot.add_line(52. + np.array([-1.4, 1.4, 1.4, -1.4, -1.4]),
              [-1.4, -1.4, 1.4, 1.4, -1.4],
              linestyle='-', linecolor='#746C66')
plot.add_line(100. + np.array([-1.4, 1.4, 1.4, -1.4, -1.4]),
              [-1.4, -1.4, 1.4, 1.4, -1.4],
              linestyle='-', linecolor='#746C66')
plot.add_line(132. + np.array([-1.4, 1.4, 1.4, -1.4, -1.4]),
              [-1.4, -1.4, 1.4, 1.4, -1.4],
              linestyle='-', linecolor='#746C66')

plot.fill_between([-15, -5, 5, 15], [-100, -100, -100, -100],
                  [-100, 100, 100, -100], fc='#F8EBC8', alpha=0.35)
plot.lines_on()
plot.markers_off()
# plot.ax.set_aspect(1.0)
plot.add_data_pointer(3., point=-100.,
                      string=r"$%.1f\,\frac{n}{s}$" % num_photoneutrons +
                      " photoneutrons emitted in $4\pi$",
                      place=(55., -85.))
plot.add_arrow(3., 25., -105., -95., fc="#7299C6")
plot.add_arrow(-3., -25., -110., -135., fc="#7299C6")
plot.add_arrow(-1., -20., -105., -101., fc="#7299C6")
plot.add_arrow(2., 10., -115., -150., fc="#7299C6")
plot.add_arrow(1., -3, -100., -75., fc="#7299C6")
plot.fill_between([-5, 5], [-5, -5], [0, 0], fc="#E3AE24", alpha=1.0,
                  leg=False)
plot.add_data_pointer(0., point=-5., place=(25., -35.),
                      string="Isocenter\n$2.4\\times 10^{4}\,\\frac{R}{hr}$")
plot.add_data_pointer(132., point=-1, place=(100., -45),
                      string=r"$>10^{2}\,\frac{R}{hr}$")
plot.add_vmeasure(0., 0., 100., offset=0.15, units='cm')
plot.add_vmeasure(0., -100., 0., offset=0.15, units='cm')
plot.add_vline(0., 52., 63., lw=0.5, color="#A7A9AC")
plot.add_arrow(0., 132., 56.5, 56.5)
plot.add_text(65., 61.5, "LS/BF3 saturated throughout")
plot.add_text(50, -115, "Concrete")
plot.add_text(50, -150, "Subfloor")
plot.add_arrow(0, 0, 100, 65, string=r"$\gamma$ beam", fc="#E3AE24")
plot.ylim(-190, 105.0)
plot.xlim(-70, 360.)
plotax = plot.ax

left, bottom = plot.fig.transFigure.inverted().transform(\
    plot.ax.transData.transform((175., -15.0)))
right, top = plot.fig.transFigure.inverted().transform(\
    plot.ax.transData.transform((350., 95.0)))
print left, bottom
print right, top
height = top - bottom
width = right - left
clinac_spectrum_ax = plot.fig.add_axes([left, bottom, width, height])
E, phi = np.loadtxt(expanduser("~") +
                    "/data/clinac_6ex_data/spectrum_bremsstrahlung.txt",
                    delimiter=',', skiprows=1, unpack=True)
plot.fill_between(E, np.zeros_like(E), 2.1E10 * np.array(phi) / 200.,
                  fc='#E3AE24', name='$\gamma$', axes=clinac_spectrum_ax)
plot.fill_between(E[E > 2.2], np.zeros_like(E[E > 2.2]),
                  2.E10 * np.array(phi[E > 2.2]) / 200.,
                  fc='#2EAFA4', name='$\gamma$', alpha=1.0,
                  axes=clinac_spectrum_ax)
plot.add_vline(2.2, 0., 4.5E11 / 200., color='#A7A9AC', axes=clinac_spectrum_ax)
plot.add_arrow(2.2, 3.5, 1.5E11 / 200., 1.5E11 / 200., '', axes=clinac_spectrum_ax)
plot.add_text(2.9, 1.65E11 / 200., r'$\left(\gamma,n\right)$ possible on $^{2}H$', axes=clinac_spectrum_ax)
plot.title("Photon spectrum at isocenter", axes=clinac_spectrum_ax)

plot.ylabel(r'Flux ($\phi$) [$\frac{\gamma}{cm^{2}s}$]', axes=clinac_spectrum_ax)
plot.xlabel(r'Energy ($E$) [$MeV$]', axes=clinac_spectrum_ax)


left, bottom = plot.fig.transFigure.inverted().transform(\
    plot.ax.transData.transform((-45., -180.)))
right, top = plot.fig.transFigure.inverted().transform(\
    plot.ax.transData.transform((-35., -110.)))
print left, bottom
print right, top
height = top - bottom
width = right - left
position = plot.fig.add_axes([left, bottom, width, height])
cb = plot.fig.colorbar(surf, cax=position, ticks=[], orientation='vertical')
position.yaxis.set_ticks_position('left')
position.yaxis.set_label_position('left')
cb.set_label(r"Photoneutron Event Density ($\rho$) [$\frac{\left(\gamma, n\right)}{cm^{3}s}$]")
left, bottom = plot.fig.transFigure.inverted().transform(\
    plot.ax.transData.transform((175., -180.0)))
right, top = plot.fig.transFigure.inverted().transform(\
    plot.ax.transData.transform((350., -75.0)))
print left, bottom
print right, top
height = top - bottom
width = right - left
spectrum_ax = plot.fig.add_axes([left, bottom, width, height])
plot.fill_between(e_0, np.zeros_like(e_0), 2.E10 * 0.000115 * n_0,
                  fc="#746C66", name='on axis', axes=spectrum_ax)
plot.fill_between(e_52, np.zeros_like(e_52), 2.E10 * 0.000115 * n_52,
                  fc="#E3AE24", name=r'$52\,cm$', axes=spectrum_ax)
plot.fill_between(e_132, np.zeros_like(e_132), 2.E10 * 0.000115 * n_132,
                  fc="#B63F97", name=r'$132\,cm$', axes=spectrum_ax, alpha=1.0)
plot.add_vline(0.5, color='#A7A9AC', axes=spectrum_ax)
plot.add_arrow(0.5, 2., 1.E1, 1.E1, '', axes=spectrum_ax)
plot.add_text(1.35, 12., "Neutron detection threshold\nfor $p_{neg}=5\,bar$", axes=spectrum_ax)

plot.ylabel(r'Flux ($\phi$) [$\frac{n}{cm^{2}s}$]', axes=spectrum_ax)
plot.xlabel(r'Energy ($E$) [$MeV$]', axes=spectrum_ax)
plot.legend(axes=spectrum_ax)
plot.title(r"Neutron spectrum from $\left(\gamma,n\right)$ interactions with floor", axes=spectrum_ax)

left, bottom = plot.fig.transFigure.inverted().transform(\
    plot.ax.transData.transform((-1.0, 0.0)))
right, top = plot.fig.transFigure.inverted().transform(\
    plot.ax.transData.transform((140.0, 50.0)))
height = top - bottom
width = right - left
n_ax = plot.fig.add_axes([left, bottom, width, height])
plot.fill_between([-1., 140.], [1.E-3, 1.E-3], [1. / 300., 1. / 300.],
                  fc='#A7A9AC', leg=False, axes=n_ax)
plot.add_text(40., 2.E-3, 'Typical Neutron Background', color="#ffffff",
              axes=n_ax)
plot = r2.plot(linecolor='#746C66', linestyle='-', addto=plot, axes=n_ax)
plot = expt.plot(linecolor='#B63F97', linestyle='-', addto=plot, axes=n_ax)

clinac_spectrum_ax.spines['right'].set_visible(False)
clinac_spectrum_ax.spines['top'].set_visible(False)
clinac_spectrum_ax.get_xaxis().tick_bottom()
clinac_spectrum_ax.get_yaxis().tick_left()
spectrum_ax.spines['right'].set_visible(False)
spectrum_ax.spines['top'].set_visible(False)
spectrum_ax.get_xaxis().tick_bottom()
spectrum_ax.get_yaxis().tick_left()
n_ax.spines['right'].set_visible(False)
n_ax.spines['top'].set_visible(False)
n_ax.get_xaxis().tick_bottom()
n_ax.get_yaxis().tick_left()

plot.xlabel(r'Distance from CLINAC Axis ($d$) [$\mathrm{cm}$]', axes=n_ax)
plot.ylabel(r'Count Rate ($\dot{c}$) [$cps$]', axes=n_ax)
plot.xlim(-1.0, 140.0, axes=n_ax)
n_ax.set_yscale("log", nonposy='clip')
plot.markers_off()
plot.lines_on()
plot.lines['experiment'].set_linewidth(0.0)
plot.lines['experiment'].set_alpha(1.0)
plot.lines['experiment'].set_markersize(6)
plot.lines['$\\frac{1}{r}$0'].set_linewidth(1.0)
plot.lines['$\\frac{1}{r}$0'].set_alpha(1.0)
plot.lines['$\\frac{1}{r}$0'].set_markersize(0)
plot.legend(axes=n_ax)


plot.ax.spines['left'].set_visible(False)
plot.ax.spines['bottom'].set_visible(False)
plot.ax.set_xticks([], [])
plot.ax.set_yticks([], [])
plot.export('../img/setup_plot', formats=['pdf', 'pgf'], sizes=['cs'],
            customsize=(14.5, 10), tight=False)

plot = ahp.ah2d()
plot.ax.set_aspect(1.0)
plot.add_circle(0, 0, 140, fc="#E3AE24", alpha=0.1)
plot.add_hmeasure(-20., 20., -30., offset=-0.01, units='cm')
plot.add_vmeasure(-20, -30., 30., units='cm')
plot.fill_between([-20, 20], [-30, -30], [30, 30], fc='#A7A9AC')
plot.add_line([0., 0.], [-35., 35.], linestyle='--')
plot.add_line([-25., 25.], [0., 0.], linestyle='--')
plot.add_data_pointer(20., point=30, place=(37.5, 40), string="Metal Table")
plot.add_data_pointer(5. / np.sqrt(2.), point=5. / np.sqrt(2.),
                      place=(37.5, 20.), string="NU or Be Target")
plot.add_circle(0, 0, 5, fc="#2EAFA4", alpha=1.0)
plot.add_data_pointer(0., point=0.,
                      place=(37.5, -20.), string="Isocenter")
plot.add_circle(0., 132., 10., fc="None", ec="#746C66")
plot.add_data_pointer(10. / np.sqrt(2.), point= 132. - 10. / np.sqrt(2.),
                      place=(35, 85), string="CTMFD Detector")
plot.add_text(55, 65, "LS/BF3 Saturated Throughout")
plot.add_hline(0., -40, -50, lw=0.5, color="#A7A9AC")
plot.add_arrow(-45, -45, 50, 132, r"$d$")
plot.add_vline(-45, 0., 45, lw=0.5, color="#A7A9AC")
plot.lines_on()
plot.markers_off()
plot.xlim(-80, 130)
plot.ylim(-60, 150)
plot.ax.spines['right'].set_visible(False)
plot.ax.spines['top'].set_visible(False)
plot.ax.spines['left'].set_visible(False)
plot.ax.spines['bottom'].set_visible(False)
plot.ax.get_xaxis().tick_bottom()
plot.ax.get_yaxis().tick_left()
plot.ax.set_xticks([], [])
plot.ax.set_yticks([], [])
plot.export('../img/clinac_setup_cartoon', formats=['pdf', 'pgf'], sizes=['cs'],
            customsize=(7.25, 7.25), tight=False)
plot.show()
