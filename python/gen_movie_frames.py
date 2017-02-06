import matplotlib
matplotlib.use("Agg")
import numpy as np
import os
import matplotlib.pyplot as plt

from interp_onto_grid import interp_onto_grid
from interp_onto_grid import plot_avg_power_2up
from c2p3_p2c3 import c2p3
from c2p3_p2c3 import p2c3

from scipy import interpolate
from scipy.spatial import Delaunay


# output grid settings
xlims = [-5, 0]
ylims = [-2.5,2.5]
zlims = [-2.5,2.5]
grid_step_size = 0.05

R_E = 6371e3

flash_lat = 31

clims = [-5, 0]
method='linear'

inp_dir = '/shared/users/asousa/WIPP/lightning_power_study/outputs/testing/lat_%d'%flash_lat
out_dir = '/shared/users/asousa/WIPP/lightning_power_study/outputs/movie_testing/'

# --------------- Latex Plot Beautification --------------------------
fig_width = 8.5 
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

# Load all files:
file_data = []
for file_ind, fname in enumerate(avail_files):
    print "loading", fname
    data_raw = np.loadtxt(os.path.join(inp_dir, fname))    
    nf = int(data_raw[0]) # num frequencies
    nv = int(data_raw[1]) # num vectors

    dlengths = np.array([int(x) for x in data_raw[2:(2 + nf*nv)]])
    nt = max(dlengths)
    dgrid = np.zeros([nf, nv, max(dlengths), 4])
    startind = 2 + nf*nv
    stopind  = 0 
    for f_ind in range(nf):
        for v_ind in range(nv):
            d_ind = f_ind*nv + v_ind
            dl = dlengths[d_ind]      
    #         print f_ind, v_ind, d_ind, dl
            stopind = startind + 4*dl
            tmp = data_raw[startind:stopind]
            curdata = np.reshape(tmp, [dl,4])
            dgrid[f_ind, v_ind, 0:dl, :] = curdata
            startind=stopind


    raw_data = dict()
    raw_data['dgrid'] = dgrid
    raw_data['nf'] = nf
    raw_data['nv'] = nv
    raw_data['nt'] = nt
    raw_data['dlengths'] = dlengths



    # Get adjacent index sets for this file:
    starting_coords = dgrid[0,:,0,0:3]

    sc_polar = np.array([c2p3(np.array(row)) for row in starting_coords])
    # Set radius to 1, round to 1-degree bins
    sc_polar[:,0] = 1
    sc_polar[:,1:]*=180./np.pi
    sc_polar = np.round(sc_polar)
    # Delaunay triangulation on lat, lon. 
    # (This is a little sloppy since we're taking cartesian distances from
    # angular values, but fuck it, it works in this case. Might break elsewhere.)
    tris = Delaunay(sc_polar[:,1:3])
    adj_inds = np.array(tris.simplices)

    raw_data['adj_inds'] = adj_inds

    file_data.append(raw_data)


print "loaded", len(file_data), "files"

tmax = max([f['nt'] for f in file_data])

print "Max time length:",tmax


xx = np.arange(xlims[0], xlims[1], grid_step_size)
yy = np.arange(ylims[0], ylims[1], grid_step_size)
zz = np.arange(zlims[0], zlims[1], grid_step_size)

nx = len(xx) 
ny = len(yy)
nz = len(zz)
print nf, nv




# Step through times (main loop):
for t in range(tmax):

    print "t:",t
    data_total = np.zeros([nx, ny, nz])     # Summed result for all frequencies (at timestep)

    # Loop over each file:
    for file in file_data:
        nf = file['nf']
        dgrid = file['dgrid']
        adj_inds =  file['adj_inds']
        dlengths = file['dlengths']

        for f_ind in range(nf):
            data = np.zeros([nx, ny, nz])
            hits = np.zeros([nx, ny, nz])
            for adj_row in adj_inds:

                # Check max length of each vector:
                cur_t = min(dlengths[adj_row])
                if t < (cur_t - 1):
                    points = dgrid[f_ind,adj_row,t:t+2,0:3]
                    vals = dgrid[f_ind,adj_row,t:t+2,3]
                    vals_flat = vals.ravel()
                    points_flat = np.vstack([points[:,:,0].ravel(), points[:,:,1].ravel(), points[:,:,2].ravel()]).T
                    
                    # Narrow down our search space to just a box around the current hull:
                    minx = min(points_flat[:,0])
                    maxx = max(points_flat[:,0])
                    miny = min(points_flat[:,1])
                    maxy = max(points_flat[:,1])
                    minz = min(points_flat[:,2])
                    maxz = max(points_flat[:,2])

                    ix = np.where((xx >= minx) & (xx <= maxx))[0]
                    iy = np.where((yy >= miny) & (yy <= maxy))[0]
                    iz = np.where((zz >= minz) & (zz <= maxz))[0]

                    px, py, pz = np.meshgrid(ix, iy, iz, indexing='ij')  # in 3d, ij gives xyz, xy gives yxz. dumb.

                    newpoints = np.vstack([xx[px.ravel()], yy[py.ravel()], zz[pz.ravel()]]).T

                    if method=='linear':
                        # If we don't specify a fill value, 'linear' mode returns NaN for anything outside
                        # the convex hull of the current point cloud (which is ideal -- we don't want any
                        # values calculated outside the (non-convex) hull.)
                        tmp_data = interpolate.griddata(points_flat, vals_flat, newpoints, method='linear', rescale=True)
                        tmp_data = tmp_data.reshape([len(ix), len(iy), len(iz)])

                        # tmp_data = np.unravel_index
                        isnans = np.isnan(tmp_data)
                        tmp_data[np.isnan(tmp_data)] = 0
            #             print np.sum(isnans), np.sum(~isnans)
                        data[px,py,pz] += tmp_data
                        hits[px,py,pz] += ~isnans
                    elif method=='mean':
                        # fill each voxel with the mean value of all corners.
                        newtri = Delaunay(points_flat)
                        hit_mask = newtri.find_simplex(newpoints) >=0

                        hit_mask = hit_mask.reshape([len(ix), len(iy), len(iz)])
                        hits[px,py,pz] += hit_mask
                        data[px[hit_mask],py[hit_mask], pz[hit_mask]] += np.mean(vals_flat)

                    # print np.max(hits), np.sum(hits > 1), len(hits)
            # Average any bins which got more than one hit at this timestep:
            # (this should just be the edges and corners)
            data[hits!=0] /= hits[hits!=0]

            data_total += data   # add the result from this frequency step to the total

        # To confirm: Total up the energy at this frame.
    total_energy = sum(sum(sum(data_total)))*pow(R_E*grid_step_size,3)
    print "Total energy:", total_energy, "Joules"

    # Plot frame:

    plot_avg_power_2up(data_total*pow(R_E*grid_step_size,3), xlims, ylims, zlims, grid_step_size, clims=clims)
    plt.savefig(os.path.join(out_dir,"frame_%d.png"%t),ldpi=300)
    plt.close('all')


