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
from scipy.integrate import nquad
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

# Volume of voxel:
def voxel_volume(points):
    tri = Delaunay(points.T, qhull_options='QJ')
    v = 0
    for row in tri.simplices:
        a = points[:,row[0]]
        b = points[:,row[1]]
        c = points[:,row[2]]
        d = points[:,row[3]]
        
        v += np.abs( np.dot(a - d, np.cross(b-d, c-d)))/6.
    return v

def rotate_latlon(raypos, itime, dlat, dlon, xf=None):
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









def interp_ray_power(ray_dir='/shared/users/asousa/WIPP/rays/2d/nightside/gcpm_kp0/',
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
                    zlims = [-2.5, 2.5],
                    step_size = 0.02,
                    clims=[1e-19, 1e-14],
                    n_sub_freqs=1,
                    frame_directory=None):

    # Constants
    Hz2Rad = 2.*np.pi
    D2R = np.pi/180.
    H_IONO_BOTTOM = 1e5
    H_IONO_TOP = 1e6

    R_E = 6371e3

    # Coordinate transform tools
    xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')

    t = np.arange(0,tmax, dt)


    # Find available rays
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
    dl = stats.mode(np.diff(lats))[0][0]

    newlons = np.array([flash_lon - dl/2., flash_lon + dl/2.])
    latgrid, longrid = np.meshgrid(lats,newlons)
    latln_pairs = zip(latgrid.ravel(), longrid.ravel())

    pairs_in_range = []


    # Adjacent frequencies to iterate over
    freqs =   [f for f in freqs if f >=f_low and f <= f_hi]
    freq_pairs = zip(freqs[0:-1],freqs[1:])


#-------------- Select points within range -----------------------------------------------------------
    for coords in latln_pairs:
        cur_d = haversine_np(flash_lon, flash_lat, coords[1], coords[0])
        if cur_d < max_dist:
            pairs_in_range.append(coords)
                
        
            
#--------------- Load and interpolate the center longitude entries ------------------------------------
    center_data = dict()
    for freq in freqs:
        print "freq: ", freq
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
            center_data[key] = curdata

#------------ Rotate center_longitude rays to new longitudes ---------------------------
    logging.info("Starting interpolation")
    ray_data = dict()
    for key in center_data.keys():
        for lon in newlons:
            newkey = (key[0], key[1], lon)
            dlon = lon - key[2] 
            d = dict()
            d['pos'] = rotate_latlon(center_data[key]['pos'],itime, 0, dlon, xf)
            d['damp']=center_data[key]['damp']
            ray_data[newkey] = d
            

# ------------- Calculate input power at each step -------------------------------------
    logging.info("Calculating input power at each cell")

    raylats =np.unique(np.array(pairs_in_range)[:,0])
    raylons =np.unique(np.array(pairs_in_range)[:,1])

    print raylats, raylons

    inp_pwrs = dict()

    lat_pairs  = zip(raylats[0:-1],raylats[1:])
    lon_pairs  = zip(raylons[0:-1],raylons[1:])
    print lon_pairs

    flash_pos_mag = [1, flash_lat, flash_lon]
    flash_pos_sm = xf.rllmag2sm(flash_pos_mag, itime)


    opts = dict()
    opts['epsabs']= 1.5e-8
    opts['epsrel']= 1.5e-8
    opts['limit']= 10

    def integrand(inlat, inlon, inw, itime, I0):
        mlt = MLT(itime, inlon, xf);
        tmp_coords = [1, inlat, inlon];
        x_sm = xf.rllmag2sm(tmp_coords, itime);

        pwr = input_power_scaling(flash_pos_sm, x_sm, inlat, inw, I0, mlt, xf);
        return pwr*(R_E + H_IONO_TOP)*D2R*(R_E + H_IONO_TOP)*np.cos(D2R*lat)*D2R

    # pwr_vec = np.zeros(len(lat_pairs))
    inp_pwrs = dict()
    for f1, f2 in freq_pairs:
        print f1, f2
        # n_freqs = np.ceil(np.abs(f2 - f1)/min_fstep)
        # f_weights = (np.arange(0,1,1.0/n_sub_freqs) + (1.0/(2.*n_sub_freqs)))
        
        # for f_weight in f_weights:
        # f_center = f_weight*f1 + (1.0-f_weight)*f2
        f_center = (f1 + f2)/2.
        for ind, (lat1, lat2) in enumerate(lat_pairs):
#             print lat1, lat2
            clat = (lat1 + lat2)/2.
            w1 = Hz2Rad*f1
            w2 = Hz2Rad*f2
            w   = Hz2Rad*(f1 + f2)/2.
            dw = np.abs(f1 - f2)*Hz2Rad
            key = (f_center, clat)
            ranges = [[lat1, lat2], lon_pairs[0]]
            # Integrate power in latitude and longitude
            integ = nquad(integrand, ranges, args=[w, itime, I0], opts=opts, full_output=False)
            pwr = integ[0]
            # Integrate in frequency (just doing this for efficiency,
            # can also do it in the nquad call)
            inp_pwrs[key] = pwr*dw


#----------- Step through and fill in the voxels (the main event) ---------------------
    logging.info("Starting interpolation")

    # output space
    xx = np.arange(xlims[0], xlims[1], step_size)
    zz = np.arange(zlims[0], zlims[1], step_size)

    nx = len(xx) 
    nz = len(zz)

    data_total = np.zeros([nx, nz])
    hits = np.zeros([nx, nz])


    interp_pos = dict()
    interp_damp= dict()

    # Interpolate between frequencies
    f_weights = (np.arange(0,1,1.0/n_sub_freqs) + (1.0/(2.*n_sub_freqs)))
    for f1, f2 in freq_pairs:
        for f_weight in f_weights:
    #         print f_weight
            tmax = 0
            for lat, lon in pairs_in_range:
                k1 = (f1, lat, lon)
                k2 = (f2, lat, lon)
                f_cur = f_weight*f1 + (1.0-f_weight)*f2
                k3 = (f_cur, lat, lon)
        #         print ray_data[k1]['nt'], ray_data[k2]['nt']
                tmax_local = min(np.shape(ray_data[k1]['pos'])[1], np.shape(ray_data[k2]['pos'])[1])
                newpos = f_weight*ray_data[k1]['pos'][:,0:tmax_local]  + (1.0 - f_weight)*ray_data[k2]['pos'][:,0:tmax_local]
                newdamp= f_weight*ray_data[k1]['damp'][0:tmax_local] + (1.0 - f_weight)*ray_data[k2]['damp'][0:tmax_local]
                interp_pos[k3] = newpos
                interp_damp[k3]=newdamp
    #             print k1, k2, k3

    lon1 = raylons[0]
    lon2 = raylons[1]
    # Step forward in time as the outer axis, so we can make movie frames:
    for t_ind in range(len(t)-1):
        # Per frequency
        data_cur = np.zeros([nx, nz])
        hits = np.zeros([nx, nz])

        print "t = ", t_ind
        for f1, f2 in freq_pairs:
            # n_freqs = np.ceil(np.abs(f2 - f1)/min_fstep)
            f_weights = (np.arange(0,1,1.0/n_sub_freqs) + (1.0/(2.*n_sub_freqs)))
    #         print f1, f2
            # Per interpolated sub-frequency between guide rays:
            for f_weight in f_weights:
                f_center = f_weight*f1 + (1.0 - f_weight)*f2
                
                # Loop over adjacent sets:
                for ind, (lat1, lat2) in enumerate(lat_pairs):
                    k0 = (f_center, lat1, lon1)
                    k1 = (f_center, lat1, lon2)
                    k2 = (f_center, lat2, lon2)
                    k3 = (f_center, lat2, lon1)
                    clat = (lat1 + lat2)/2.

                    tmax_local = min(len(interp_damp[k0]),len(interp_damp[k1]),
                                     len(interp_damp[k2]), len(interp_damp[k3]))
                    if t_ind < tmax_local - 1:
                        points = np.hstack([interp_pos[k0][:,t_ind:t_ind+2], 
                                            interp_pos[k1][:,t_ind:t_ind+2], 
                                            interp_pos[k2][:,t_ind:t_ind+2],
                                            interp_pos[k3][:,t_ind:t_ind+2]])
                        damps  = np.hstack([interp_damp[k0][t_ind:t_ind+2], 
                                            interp_damp[k1][t_ind:t_ind+2], 
                                            interp_damp[k2][t_ind:t_ind+2],
                                            interp_damp[k3][t_ind:t_ind+2]])
                        pwr = inp_pwrs[((f1 + f2)/2., clat)]/n_sub_freqs

    #                     print t_ind, np.shape(points)
                        avgy = np.mean(points[1,:])
                        # This block for Cartesian-gridded output space:
                        minx = min(points[0,:])
                        maxx = max(points[0,:])
                        minz = min(points[2,:])
                        maxz = max(points[2,:])
                        ix = np.where((xx >= minx) & (xx <= maxx))[0]
                        iz = np.where((zz >= minz) & (zz <= maxz))[0]
                        px, pz = np.meshgrid(ix, iz, indexing='ij')  # in 3d, ij gives xyz, xy gives yxz. dumb.
                        newpoints = np.vstack([xx[px.ravel()],np.ones_like(px.ravel())*avgy, zz[pz.ravel()]])
                        damping_avg = np.mean(damps)
                        tri = Delaunay(points.T,qhull_options='QJ')
                        mask = (tri.find_simplex(newpoints.T) >= 0)*1.0
                        mask = mask.reshape([len(ix), len(iz)])
                        total_cells = np.sum(mask)

                        if (total_cells > 0):
                            voxel_vol = voxel_volume(points) # Volume in R_e           
                            data_cur[px,pz] += damping_avg*mask*pwr/voxel_vol/pow(R_E,3)  # better volume estimate   
        
        # Plot individual timesteps (for movies)
        if frame_directory is not None:
            fig, ax = plot_xz(data_cur, xlims,zlims, step_size, clims)
            fig.suptitle('T=%0.2f s, %0.1f J'%(t_ind*dt, np.sum(data_cur)*pow(step_size*R_E,3)))
            fig.savefig(os.path.join(frame_directory,'frame%d.png'%t_ind),ldpi=300)
            plt.close('all')

        data_total += data_cur
            
    logging.info("finished with interpolation")
    return data_total
    # plot_xz(data_total, xlims, zlims, step_size)

def plot_xz(data, xlims, zlims, step_size, clims=None):
    # --------------- Latex Plot Beautification --------------------------
    fig_width = 6 
    fig_height = 6
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

    # Constants
    Hz2Rad = 2.*np.pi
    D2R = np.pi/180.
    H_IONO_BOTTOM = 1e5
    H_IONO_TOP = 1e6

    R_E = 6371e3


    xx = np.arange(xlims[0], xlims[1], step_size)
    zz = np.arange(zlims[0], zlims[1], step_size)

    nx = len(xx) 
    nz = len(zz)

    logdata = np.log10(data)
    logdata[np.isinf(logdata)] = -100

    maxlog = np.max([logdata])
    
    # Show about 5 orders of magnitude
    if clims is None:
        clims = [maxlog - 5, maxlog]
    
    fig, ax = plt.subplots(1,1)
    # Plot the earth
    earth = plt.Circle((0,0),1,color='0.5',alpha=1, zorder=100)
    iono  = plt.Circle((0,0),(R_E + H_IONO_TOP)/R_E, color='w',alpha=0.8, zorder=99)
    ax.add_patch(earth)   
    ax.add_patch(iono)
    
    p0 = ax.pcolorfast(xx, zz, logdata.T, vmin=clims[0], vmax=clims[1])
    ax.set_aspect('equal')
    ax.set_xlim([xx[0],xx[-1]])
    ax.set_ylim([zz[0],zz[-1]])

    fig.tight_layout()
    
    fig.subplots_adjust(right=0.83)
    cax = fig.add_axes([0.85,0.14, 0.02, 0.75])

    cb = plt.colorbar(p0, cax=cax)
    cb.set_label('avg wave power density')
    cticks = np.arange(clims[0],clims[1] + 1)
    cb.set_ticks(cticks)
    cticklabels = ['$10^{%d}$'%k for k in cticks]
    cb.set_ticklabels(cticklabels)

    return fig, ax









if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO,
                        format='[%(levelname)s] %(message)s')  


    xlims = [-8, 0]
    zlims = [-4, 4]

    datagrid = interp_ray_power(ray_dir='/shared/users/asousa/WIPP/rays/2d/nightside/gcpm_kp0',
                                tmax = 15,
                                flash_lat=40,
                                flash_lon=76,
                                dt=0.1,
                                f_low=200,
                                f_hi=30000,
                                max_dist=1500,
                                xlims = xlims,
                                zlims = zlims,
                                n_sub_freqs=50,
                                clims=[-19, -13],
                                frame_directory='/shared/users/asousa/WIPP/lightning_power_study/outputs/testing'
                                )
    # datagrid = interp_ray_power()


    fig, ax = plot_xz(datagrid, 
                xlims=xlims,
                zlims=zlims,
                step_size=0.02)

    fig.savefig("test_figure.png",ldpi=300)
    
            

 