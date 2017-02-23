import matplotlib
matplotlib.use('agg')
import numpy as np
import pandas as pd
import pickle
from scipy import interpolate
import matplotlib.pyplot as plt
import os
import itertools
import random
import os
import time
import datetime as datetime
from raytracer_utils import read_rayfile, read_damp
from scipy.spatial import Delaunay
from scipy import stats
import xflib
from graf_iono_absorp import total_input_power, MLT, input_power_scaling
import logging





def haversine_np(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    All args must be of equal length.    

    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6371 * c
    return km

def voxel_volume(points):
    '''
    Compute the volume of a voxel given by an array of points
    '''
    tri = Delaunay(points.T)
    v = 0
    for row in tri.simplices:
        a = points[:,row[0]]
        b = points[:,row[1]]
        c = points[:,row[2]]
        d = points[:,row[3]]
        
        v += np.abs( np.dot(a - d, np.cross(b-d, c-d)))/6.
    return v


def rotate_latlon(raypos, itime, dlat, dlon, xf=None):
    ''' 1. Convert a ray from SM to rllmag
        2. Rotate by dlat and dlon
        3. Convert back to SM and return a copy
    '''
    if xf is None:
        xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')
        
    newpos = np.zeros_like(raypos)
    for ind in range(np.shape(raypos)[1]):
#         print ind
        tmp = xf.sm2rllmag(raypos[:,ind], itime)
        tmp[1] += dlat
        tmp[2] += dlon
        newpos[:,ind] = xf.rllmag2sm(tmp, itime)
    
    return newpos


def interp_ray_power(ray_dir='/shared/users/asousa/WIPP/rays/2d/ngo_dipole/',
                    flash_lat=40,
                    flash_lon=76,
                    max_dist=1500, 
                    tmax=10,
                    dt=0.1,
                    I0=-10000,
                    f_low=200, f_hi=30000,
                    min_fstep = 5000.0,
                    itime = datetime.datetime(2010,1,1,0,0,0),
                    xlims = [-5, 0],
                    ylims = [-2.5, 2.5],
                    zlims = [-2.5, 2.5],
                    step_size = 0.02,
                    frame_directory=None):

    n_freqs = 1  # Move this somewhere else, and make it a function of the guide frequencies

    # Constants
    Hz2Rad = 2.*np.pi
    D2R = np.pi/180.
    H_IONO = 1e5
    R_E = 6371e3

    # Coordinate transform tools
    xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')

    # New time axis
    t = np.arange(0,tmax, dt)

    # Find available files
    d = os.listdir(ray_dir)
    freqs = sorted([int(f[2:]) for f in d if f.startswith('f_')])
    d = os.listdir(os.path.join(ray_dir, 'f_%d'%freqs[0]))
    lons = sorted([float(f[4:]) for f in d if f.startswith('lon_')])
    d = os.listdir(os.path.join(ray_dir, 'f_%d'%freqs[0], 'lon_%d'%lons[0]))
    lats = sorted([float(s.split('_')[2]) for s in d if s.startswith('ray_')])

    # Closest available longitude to the flash:
    center_lon = lons[np.argmin(np.abs(np.array(lons) - flash_lon))]
    center_lat = lats[np.argmin(np.abs(np.array(lats) - flash_lat))]

    # Latitude spacing:
    # dl = stats.mode(np.diff(lats))[0][0]
    dl = 1

    newlons = np.arange(center_lon - 2*np.round(max_dist/111.0), center_lon + 2*np.round(max_dist/111.0) + dl, dl, dtype=float)
    latgrid, longrid = np.meshgrid(lats,newlons)
    latln_pairs = zip(latgrid.ravel(), longrid.ravel())

    pairs_in_range = []

    # Adjacent frequencies to iterate over
    freqs =   [f for f in freqs if f >=f_low and f <= f_hi]
    freq_pairs = zip(freqs[0:-1],freqs[1:])



#-------------- Select points within range -----------------------------------------------------------
    # Prune out some points the further out we go:
    logging.info("selecting guide rays")
    for coords in latln_pairs:
        if coords[0]%1==0:
            cur_d = haversine_np(flash_lon, flash_lat, coords[1], coords[0])
            if cur_d < max_dist:

                if (cur_d < 300):
                    pairs_in_range.append(coords)
                elif (cur_d < 600 and (center_lat - coords[0])%2==1 and (center_lon - coords[1])%2==1):
                    pairs_in_range.append(coords)
                elif (cur_d < 1200 and (center_lat - coords[0])%4==1 and (center_lon - coords[1])%4==1):                
                    pairs_in_range.append(coords)
                elif (cur_d < 2000 and (center_lat - coords[0])%8==1 and (center_lon - coords[1])%8==1):                
                    pairs_in_range.append(coords)

            
#--------------- Load and interpolate the center longitude entries ------------------------------------
    logging.info("Loading center-longitude rays")

    ray_data = dict()
    for freq in freqs:
        for lat in np.unique([x[0] for x in pairs_in_range]):
            lon = center_lon
            filename = os.path.join(ray_dir,'f_%d'%freq,'lon_%d'%lon,'ray_%d_%d_%d.ray'%(freq,lat,lon))
    #         print filename
            rf = read_rayfile(filename)[0]
            
            filename = os.path.join(ray_dir,'f_%d'%freq,'lon_%d'%lon,'damp_%d_%d_%d.ray'%(freq,lat,lon))
            df = read_damp(filename)[0]
            
            t_cur = t[t <= rf['time'].iloc[-1]]
            
            # Interpolate onto our new time axis:
            x = interpolate.interp1d(rf['time'],rf['pos']['x']).__call__(t_cur)/R_E
            y = interpolate.interp1d(rf['time'],rf['pos']['y']).__call__(t_cur)/R_E
            z = interpolate.interp1d(rf['time'],rf['pos']['z']).__call__(t_cur)/R_E
            d = interpolate.interp1d(df['time'],df['damping']).__call__(t_cur)

            # Stash it somewhere:
            key = (freq, lat, lon)
            curdata = dict()
            curdata['pos'] = np.vstack([x,y,z])
            curdata['damp']= d
            curdata['nt'] = len(t_cur)
            ray_data[key] = curdata




# ------------- Get adjacent triangles -------------------------------------------------------------------
    logging.info("Delaunay triangualtion")

    tris = Delaunay(pairs_in_range)
    adj_inds = np.array(tris.simplices)
    starting_coords = np.array(pairs_in_range)

    logging.info("%d simplexes"%np.shape(adj_inds)[0])


    # # Loop through frequency pairs:
    # f1 = 3920
    # f2 = 4580


    # print min(starting_coords[:,0]), max(starting_coords[:,0])

    grid_spacing = 0.25 # deg (Just for the power calculation)

    latrange = [min(starting_coords[:,0]), max(starting_coords[:,0])]
    lonrange = [min(starting_coords[:,1]), max(starting_coords[:,1])]
    gridlats = np.arange(latrange[0]-grid_spacing/2, latrange[1] + grid_spacing/2, grid_spacing)
    gridlons = np.arange(lonrange[0]-grid_spacing/2, lonrange[1] + grid_spacing/2, grid_spacing)

    clats, clons = np.meshgrid(gridlats, gridlons)

    mask = tris.find_simplex(zip(clats.ravel(), clons.ravel()))
    # print np.shape(mask)
    mask = np.reshape(mask, [len(gridlons), len(gridlats)])
    # print np.shape(mask)
    # print np.shape(clats)


    flash_pos_mag = [1, flash_lat, flash_lon]
    flash_pos_sm = xf.rllmag2sm(flash_pos_mag, itime)

    logging.info("Calculating input power at each simplex")

    inp_pwrs = dict()

    for f1, f2 in freq_pairs:
        n_freqs = np.ceil(np.abs(f2 - f1)/min_fstep)
        f_weights = (np.arange(0,1,1.0/n_freqs) + (1.0/(2.*n_freqs)))
        # logging.info("f1: %g, f2: %g"%(f1,f2))

        for f_weight in f_weights:
            f_center = f_weight*f1 + (1.0-f_weight)*f2
            logging.info("freq: %g"%f_center)
            pwr_grid = np.zeros([len(gridlons), len(gridlats)])
            # Calculate power within enclosed triangles:
            for lat_ind, lat in enumerate(gridlats):
                for lon_ind, lon in enumerate(gridlons):
                    if mask[lon_ind, lat_ind] >=0:
                            clat = lat + grid_spacing/2.
                            clon = lon + grid_spacing/2.
                            w = f_center*Hz2Rad
                            # w = np.abs(f1 + f2)*Hz2Rad/2.;
                            dw = np.abs(f1 - f2)*Hz2Rad/n_freqs
                            mlt = MLT(itime, clon, xf);
                            tmp_coords = [1 + H_IONO/R_E, clat, clon];
                            x_sm = xf.rllmag2sm(tmp_coords, itime);
                            pwr = input_power_scaling(flash_pos_sm, x_sm, lat, w, I0, mlt, xf);
                            dist_lat = (R_E + H_IONO)*grid_spacing*D2R;
                            dist_lon = (R_E + H_IONO)*grid_spacing*np.cos(D2R*clat)*D2R;
                            pwr_grid[lon_ind, lat_ind] = pwr*dw*dist_lat*dist_lon
                                    

            tri_inds = np.unique(mask[mask >= 0])
            # inp_pwrs = np.zeros(len(tri_inds))

            for val in tri_inds:
                inp_pwrs[(f_center, val)] = np.sum(pwr_grid[mask==val])
                # inp_pwrs[val] = np.sum(pwr_grid[mask==val])
                        
    # Final outputs of this section:
    # inp_pwrs ~ total energy (joules) within each triangle
    # adj_inds ~ indices of corner vectors


#------------ Rotate center_longitude rays to new longitudes ---------------------------
    logging.info("rotating center longitude rays to offset longitudes")
    for freq in freqs:
        for lat, lon in pairs_in_range:
            key = (freq, lat, lon)
            if not key in ray_data:
                centerkey = (freq, lat, center_lon)
                centerray = ray_data[centerkey]
                dlon = lon - center_lon
                d = dict()
                
                d['pos'] = rotate_latlon(centerray['pos'],itime, 0, dlon, xf)
                d['damp']= centerray['damp']
                d['nt'] = centerray['nt']
                ray_data[key] = d
                


#----------- Step through and fill in the voxels (the main event) ---------------------
    logging.info("Starting interpolation")
    # output space
    xx = np.arange(xlims[0], xlims[1], step_size)
    yy = np.arange(ylims[0], ylims[1], step_size)
    zz = np.arange(zlims[0], zlims[1], step_size)

    nx = len(xx) 
    ny = len(yy)
    nz = len(zz)

    data_total = np.zeros([nx, ny, nz])
    hits = np.zeros([nx, ny, nz])


    # Interpolate over frequencies:
    # freq_step = (f2 - f1)/n_freqs
    
    interp_pos = dict()
    interp_damp= dict()


    # for f_weight in np.linspace(0,1, n_freqs):
    # f_weight is the center of each of n_freqs subsections between 0 and 1
    for f1, f2 in freq_pairs:
        n_freqs = np.ceil(np.abs(f2 - f1)/min_fstep)
        f_weights = (np.arange(0,1,1.0/n_freqs) + (1.0/(2.*n_freqs)))
        for f_weight in f_weights:
            logging.info("Frequency weight: %g"%f_weight)
            tmax = 0
            # Interpolate in frequency between guide rays
            for lat, lon in pairs_in_range:
                k1 = (f1, lat, lon)     # Lower guide ray
                k2 = (f2, lat, lon)     # Upper guide ray
                f_cur = f_weight*f1 + (1.0-f_weight)*f2
                k3 = (f_cur, lat, lon)  # Interpolated ray
                
                tmax_local = min(ray_data[k1]['nt'], ray_data[k2]['nt'])
                newpos = f_weight*ray_data[k1]['pos'][:,0:tmax_local]  + (1.0 - f_weight)*ray_data[k2]['pos'][:,0:tmax_local]
                newdamp= f_weight*ray_data[k1]['damp'][0:tmax_local] + (1.0 - f_weight)*ray_data[k2]['damp'][0:tmax_local]
                interp_pos[k3] = newpos
                interp_damp[k3]=newdamp

#     Step through times:
    for t_ind in range(len(t)-1):
        logging.info("t= %g sec"%(t_ind*dt)) 
        data_cur = np.zeros([nx, ny, nz])
        for f1, f2 in freq_pairs:
            n_freqs = np.ceil(np.abs(f2 - f1)/min_fstep)
            f_weights = (np.arange(0,1,1.0/n_freqs) + (1.0/(2.*n_freqs)))

            for f_weight in f_weights:
                f_center = f_weight*f1 + (1.0 - f_weight)*f2
                logging.info("Frequency: %g"%f_center)

                # Loop over adjacent sets:
                for ind, adj_row in enumerate(adj_inds):
                    k0 = (f_center, pairs_in_range[adj_row[0]][0], pairs_in_range[adj_row[0]][1])
                    k1 = (f_center, pairs_in_range[adj_row[1]][0], pairs_in_range[adj_row[1]][1])
                    k2 = (f_center, pairs_in_range[adj_row[2]][0], pairs_in_range[adj_row[2]][1])

                    tmax_local = min(len(interp_damp[k0]),len(interp_damp[k1]),len(interp_damp[k2]))
                    if t_ind < tmax_local - 1:
                        points = np.hstack([interp_pos[k0][:,t_ind:t_ind+2], 
                                            interp_pos[k1][:,t_ind:t_ind+2], 
                                            interp_pos[k2][:,t_ind:t_ind+2]])
                        damps  = np.hstack([interp_damp[k0][t_ind:t_ind+2], 
                                            interp_damp[k1][t_ind:t_ind+2], 
                                            interp_damp[k2][t_ind:t_ind+2]])

                        pwr = inp_pwrs[(f_center, ind)]

                        minx = min(points[0,:])
                        maxx = max(points[0,:])
                        miny = min(points[1,:])
                        maxy = max(points[1,:])
                        minz = min(points[2,:])
                        maxz = max(points[2,:])

                        ix = np.where((xx >= minx) & (xx <= maxx))[0]
                        iy = np.where((yy >= miny) & (yy <= maxy))[0]
                        iz = np.where((zz >= minz) & (zz <= maxz))[0]

                        px, py, pz = np.meshgrid(ix, iy, iz, indexing='ij')  # in 3d, ij gives xyz, xy gives yxz. dumb.
                        newpoints = np.vstack([xx[px.ravel()], yy[py.ravel()], zz[pz.ravel()]]).T

                        # Average the damping and just fill in the voxel:
                        damping_avg = np.mean(damps)
                        tri = Delaunay(points.T)
                        mask = (tri.find_simplex(newpoints) >= 0)*1.0
                        mask = mask.reshape([len(ix), len(iy), len(iz)])
                        total_cells = np.sum(mask)
                        if (total_cells > 0):
            #                voxel_vol = voxel_volume(points) # Volume in R_e           
            #                 data_cur[px,py,pz] += damping_avg*mask*pwr/voxel_vol/pow(R_E,3)  # better volume estimate   
                            data_cur[px,py,pz] += damping_avg*mask*pwr/total_cells/pow(R_E*step_size,3)  # assures constant energy

                    # Average any bins which got more than one hit at this timestep:
                    # (this should just be the edges and corners)
            #         data_cur[hits!=0] /= hits[hits!=0]

        # Plot individual timesteps (for movies)
        if frame_directory is not None:
            fig, ax = plot_avg_pwr(data_cur, xlims,ylims,zlims,step_size)
            total_energy = sum(sum(sum(data_cur)))*pow(R_E*step_size,3)
            fig.suptitle('T=%0.f s E= %0.1f J'%(t_ind*dt, total_energy))
            fig.savefig(os.path.join(frame_directory,'frame%d.png'%t_ind),ldpi=300)
            plt.close('all')


        # logging.info("t ",t_ind, "total energy ", np.sum((data_cur))*pow(R_E*step_size,3), "multi hits: ", np.sum(hits > 1))
        data_total += data_cur

        logging.info("energy: %g"%(np.sum((data_cur))*pow(R_E*step_size,3)))
    logging.info("finished with interpolation")
    return data_total



def plot_avg_pwr(data, xlims, ylims, zlims, step_size):
    ''' 3-up plot.
    '''
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
    
    R_E = 6371e3
    H_IONO = 1000e3

    logging.info('Plotting')
    xx = np.arange(xlims[0], xlims[1], step_size)
    yy = np.arange(ylims[0], ylims[1], step_size)
    zz = np.arange(zlims[0], zlims[1], step_size)

    nx = len(xx) 
    ny = len(yy)
    nz = len(zz)

    xz_sum = np.log10((np.sum(data, axis=1).T)/len(yy))
    xy_sum = np.log10((np.sum(data, axis=2).T)/len(zz))
    yz_sum = np.log10((np.sum(data, axis=0).T)/len(xx))

    maxlog = np.max([np.max(xz_sum), np.max(xy_sum), np.max(yz_sum)])
    
    # Show about 3 orders of magnitude
    clims = [maxlog - 3, maxlog]

    xy_sum[np.isinf(xy_sum)] = -100
    xz_sum[np.isinf(xz_sum)] = -100
#     yz_sum[np.isinf(yz_sum)] = -100
    flatblack = np.ones_like(yz_sum)*-100

    fig, ax = plt.subplots(1,3)
    # Plot the earth
    for i in [0, 1, 2]:
        earth = plt.Circle((0,0),1,color='0.5',alpha=1, zorder=100)
        iono  = plt.Circle((0,0),(R_E + H_IONO)/R_E, color='w',alpha=0.8, zorder=99)
        ax[i].add_patch(earth)   
        ax[i].add_patch(iono)
    
    p0 = ax[0].pcolorfast(xx, yy, xy_sum, vmin=clims[0], vmax=clims[1])
    p1 = ax[1].pcolorfast(xx, zz, xz_sum, vmin=clims[0], vmax=clims[1])
    p2 = ax[2].pcolorfast(yy, zz, yz_sum, vmin=clims[0], vmax=clims[1], zorder=101)
    p3 = ax[2].pcolorfast(yy, zz, flatblack, vmin=clims[0], vmax=clims[1], zorder=98)
    
    ax[0].set_aspect('equal')
    ax[1].set_aspect('equal')
    ax[2].set_aspect('equal')

    ax[0].set_title('XY')
    ax[1].set_title('XZ')
    ax[2].set_title('YZ')
    
    ax[0].set_xlim([xx[0],xx[-1]])
    ax[0].set_ylim([yy[0],yy[-1]])
    ax[1].set_xlim([xx[0],xx[-1]])
    ax[1].set_ylim([zz[0],zz[-1]])
    ax[2].set_xlim([yy[0],yy[-1]])
    ax[2].set_ylim([zz[0],zz[-1]])
    
        
    fig.tight_layout()
    
    
    # cax = fig.colorbar(p, ax=ax.ravel().tolist())
    fig.subplots_adjust(right=0.84)
    cax = fig.add_axes([0.85,0.19, 0.01, 0.629])

    cb = plt.colorbar(p2, cax=cax)
    cb.set_label('avg wave power density')
    cticks = np.arange(clims[0],clims[1] + 1)
    cb.set_ticks(cticks)
    cticklabels = ['$10^{%d}$'%k for k in cticks]
    cb.set_ticklabels(cticklabels)
    
    return fig, ax




if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO,
                        format='[%(levelname)s] %(message)s')  

    datagrid = interp_ray_power(ray_dir='/shared/users/asousa/WIPP/rays/2d/ngo_dipole',
                                tmax = 15,
                                flash_lat=45,
                                flash_lon=77,
                                dt=0.05,
                                frame_directory='/shared/users/asousa/WIPP/lightning_power_study/outputs/movie_testing7/')
    # datagrid = interp_ray_power()


    fig, ax = plot_avg_pwr(datagrid, 
                xlims=[-5,0],
                ylims=[-2.5,2.5],
                zlims=[-2.5,2.5],
                step_size=0.02)

    fig.savefig("test_figure.png",ldpi=300)
    