import matplotlib
matplotlib.use("Agg")
import numpy as np
import os
import matplotlib.pyplot as plt

from interp_onto_grid import interp_onto_grid
from interp_onto_grid import plot_avg_power_3up



# output grid settings
xlims = [-5, 0]
ylims = [-2.5,2.5]
zlims = [-2.5,2.5]
grid_step_size = 0.02

flash_lat = 50

inp_dir = '/shared/users/asousa/WIPP/lightning_power_study/outputs/nightside/ngo_dipole/lat_%d'%flash_lat
# out_dir = '/shared/users/asousa/WIPP/lightning_power_study/outputs/nightside/ngo_dipole/lat_%d'%flash_lat
out_dir = '/shared/users/asousa/WIPP/lightning_power_study/outputs/nightside/test_figs'



# --------------- Latex Plot Beautification --------------------------
fig_width = 12 
fig_height = 4
fig_size =  [fig_width+1,fig_height+1]
params = {'backend': 'ps',
          'axes.labelsize': 14,
          'font.size': 14,
          'legend.fontsize': 14,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14,
          'text.usetex': False,
          'figure.figsize': fig_size}
plt.rcParams.update(params)
# --------------- Latex Plot Beautification --------------------------




d = os.listdir(inp_dir)
avail_files = [x for x in d if x.startswith('power_vectors') and x.endswith('.dat')]


data_total = None

for file in avail_files:
    print "file: ", file
    curfile = os.path.join(inp_dir, file)
    data = interp_onto_grid(curfile, xlims, ylims, zlims, grid_step_size, method='linear')
    plot_avg_power_3up(data, xlims,ylims,zlims,grid_step_size)
    figfile = os.path.join(out_dir, "figure_" + file + ".png")
    plt.savefig(figfile, ldpi=300)
    plt.close('all')

    # dumpfile = os.path.join(out_dir, "dump_" + file)
    # np.save(dumpfile, data)

    if data_total is None:
        data_total = data
    else:
        data_total += data


plot_avg_power_3up(data_total, xlims, ylims, zlims, grid_step_size, clims=[-12,-5])
plt.savefig(os.path.join(out_dir,"figure_TOTAL.png"),ldpi=300)
plt.close('all')

dumpfile = os.path.join(out_dir, "dump_TOTAL")
np.save(dumpfile, data_total)
