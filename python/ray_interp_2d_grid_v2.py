# import matplotlib
# matplotlib.use('agg')
import numpy as np
import pandas as pd
import pickle
import gzip
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
from graf_iono_absorp import total_input_power, lon2MLT, MLT2lon, input_power_scaling
import logging
import math


# ----------------------------- Methods -------------------------------
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

def rotate_latlon(raypos, itime, dlat, dlon, xf=None, output = 'SM'):
    if xf is None:
        xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')
        
    newpos = np.zeros_like(raypos)
    for ind in range(np.shape(raypos)[1]):
#         print ind
        tmp = xf.sm2rllmag(raypos[:,ind], itime)
        tmp[1] += dlat
        tmp[2] += dlon
        if output=='SM':
            newpos[:,ind] = xf.rllmag2sm(tmp, itime)
        elif output=='MAG':
            newpos[:,ind] = xf.s2c(tmp)
    return newpos

def flatten_longitude_variation(raypos, itime, xf=None):
    if xf is None:
        xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')
        
    newpos = np.zeros_like(raypos)

    tmp = xf.sm2rllmag(raypos[:,0], itime)
    start_lon = tmp[2]

    for ind in range(np.shape(raypos)[1]):
#         print ind
        tmp = xf.sm2rllmag(raypos[:,ind], itime)
        # tmp[1] += dlat
        tmp[2] = start_lon
        newpos[:,ind] = xf.rllmag2sm(tmp, itime)
    
    return newpos


def voxel_vol_nd(points):
    '''
    volume of a polygon in n-dimensional space. Rad.
    '''
    n, m = np.shape(points)
    tri = Delaunay(points.T, qhull_options='QJ')
    v = 0
    for row in tri.simplices:
        mat = points[:, row[1:]].T - points[:, row[0]].T
        v += np.abs(np.linalg.det(mat)/math.factorial(n))
    return v

def interp_ray_power(ray_dir=None,
                    power_dir = None,
                    flash_lat=40,
                    mlt = 0,
                    max_dist=1500, 
                    tmax=10,
                    dt=0.1,
                    I0=-10000,
                    f_low=200, f_hi=30000,
                    min_fstep = 5000.0,
                    itime = datetime.datetime(2010,1,1,0,0,0),
                    step_size = 0.02,
                    clims=[1e-19, 1e-14],
                    n_sub_freqs=1,
                    dlon = 1,
                    Llims = [0, 5],
                    NL = 100,
                    frame_directory=None):
    ''' This version is a re-do of the 2d grid interpolation, using all the 
        changes I made for the L-shell version. Intending this to be the version
        used for movie frames in the defense slides.'''

    # Constants
    Hz2Rad = 2.*np.pi
    D2R = np.pi/180.
    H_IONO_BOTTOM = 1e5
    H_IONO_TOP = 1e6
    R_E = 6371e3






    # Coordinate transform tools
    xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')
    # time axis
    t = np.arange(0,tmax, dt)
    itime = datetime.datetime(2010,1,1,0,0,0)

    # Select flash longitude corresponding to the input MLT:
    flash_lon = MLT2lon(itime,mlt,xf)

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
    dlat = stats.mode(np.diff(lats))[0][0]

    newlons = np.array([flash_lon - dlon/2., flash_lon + dlon/2.])
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
        # logging.info("Loading freq: %d"%freq)
        print "Loading freq: %d"%freq
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
            d = interpolate.interp1d(df['time'],df['damping'], bounds_error=False, fill_value=0).__call__(t_cur)
            
            # Stash it somewhere:
            key = (freq, lat, lon)
            curdata = dict()

            # Flatten out any longitude variation, just to be sure:
            curdata['pos'] = flatten_longitude_variation(np.vstack([x,y,z]), itime, xf=xf)
            # curdata['pos'] = np.vstack([x,y,z])
            curdata['damp']= d
            curdata['nt'] = len(t_cur)
            center_data[key] = curdata

#------------ Rotate center_longitude rays to new longitudes ---------------------------
    # logging.info("Rotating to new longitudes")
    print "Rotating to new longitudes"

    ray_data = dict()
    for key in center_data.keys():
        for lon in newlons:
            newkey = (key[0], key[1], lon)
            dlon = lon - key[2] 
            d = dict()
            d['pos'] = rotate_latlon(center_data[key]['pos'],itime, 0, dlon, xf, output='MAG')
            d['damp']=center_data[key]['damp']
            ray_data[newkey] = d


# ------------- Calculate input power at each step -------------------------------------
    # logging.info("Calculating input power at each cell")

    raylats =np.unique(np.array(pairs_in_range)[:,0])
    raylons =np.unique(np.array(pairs_in_range)[:,1])

    # print raylats, raylons


    lat_pairs  = zip(raylats[0:-1],raylats[1:])
    # lon_pairs  = zip(flash_lons[0:-1],flash_lons[1:])
    # # print lon_pairs

    flash_pos_mag = [1, flash_lat, flash_lon]
    flash_pos_sm = xf.rllmag2sm(flash_pos_mag, itime)

    try:
        print "Loading previous powers"
        with gzip.open(os.path.join(power_dir,'input_energy_%d_%d.pklz'%(flash_lat, mlt)),'r') as file:
            pwr_db = pickle.load(file)
    except:
        print "Failed to load power file"


#----------- Step through and fill in the voxels (the main event) ---------------------
    # logging.info("Starting interpolation")
    print "Starting interpolation"
    # output space

    # Set up output grid:

    xx = np.linspace(Llims[0], Llims[1], NL)
    zz = np.linspace(-Llims[1]/2., Llims[1]/2., NL)


    nt = len(t)
    n_freq_pairs = len(freq_pairs)
    # data_total = np.zeros([NL, NL, n_freq_pairs, nt])
    data_total = np.zeros([NL, NL, nt])

    lon1 = raylons[0]
    lon2 = raylons[1]

    for t_ind in np.arange(nt - 1):
        # Per frequency
        print "t = ", t_ind
        data_cur = np.zeros([NL, NL])
        for freq_ind, (f1, f2) in enumerate(freq_pairs):

            # Fine-scale frequencies to use
            ff = np.arange(0, n_sub_freqs, 1)
            nf = len(ff)
            dc = np.zeros([NL, NL, nf])

            # Loop over each pair of adjacent latitudes
            for lat1, lat2 in lat_pairs:
                # dc = np.zeros([NL, NL, nf])
                # --- Calculate Delaunay triangulation given the current timestep:
                k0 = (f1, lat1, lon1)
                k1 = (f1, lat2, lon1)
                k2 = (f2, lat1, lon1)
                k3 = (f2, lat2, lon1)
                k4 = (f1, lat1, lon2)
                k5 = (f1, lat2, lon2)
                k6 = (f2, lat1, lon2)
                k7 = (f2, lat2, lon2)
                clat = (lat1 + lat2)/2.
                f_center = (f1 + f2)/2.

                tmax_local = min(np.shape(ray_data[k0]['pos'])[1], np.shape(ray_data[k1]['pos'])[1],
                                 np.shape(ray_data[k2]['pos'])[1], np.shape(ray_data[k3]['pos'])[1],
                                 np.shape(ray_data[k4]['pos'])[1], np.shape(ray_data[k5]['pos'])[1],
                                 np.shape(ray_data[k6]['pos'])[1], np.shape(ray_data[k7]['pos'])[1])
                # Don't do it if we're past the end of the shortest ray
                if (t_ind < tmax_local - 1):

                    points_4d = np.hstack([np.vstack([ray_data[k0]['pos'][:,t_ind:t_ind+2],np.zeros([1,2])]),
                                           np.vstack([ray_data[k1]['pos'][:,t_ind:t_ind+2],np.zeros([1,2])]),
                                           np.vstack([ray_data[k2]['pos'][:,t_ind:t_ind+2],np.ones([1,2])*nf]),
                                           np.vstack([ray_data[k3]['pos'][:,t_ind:t_ind+2],np.ones([1,2])*nf]),
                                           np.vstack([ray_data[k4]['pos'][:,t_ind:t_ind+2],np.zeros([1,2])]),
                                           np.vstack([ray_data[k5]['pos'][:,t_ind:t_ind+2],np.zeros([1,2])]),
                                           np.vstack([ray_data[k6]['pos'][:,t_ind:t_ind+2],np.ones([1,2])*nf]),
                                           np.vstack([ray_data[k7]['pos'][:,t_ind:t_ind+2],np.ones([1,2])*nf])])
                
                    voxel_vol = voxel_vol_nd(points_4d)*pow(R_E,3.)

                    damps_2d = np.hstack([ray_data[k0]['damp'][t_ind:t_ind+2],
                                          ray_data[k1]['damp'][t_ind:t_ind+2],
                                          ray_data[k2]['damp'][t_ind:t_ind+2],
                                          ray_data[k3]['damp'][t_ind:t_ind+2]])
                    damping_avg = np.mean(damps_2d)
                    
                    # pwr = inp_pwrs[(f_center, clat)][0]  # Avg pwrs per frequency bin (Take center longitude)
                    pwr_freq_ind = np.argmin(np.abs(pwr_db['cfreqs'] - f_center))
                    pwr_lat_ind  = np.argmin(np.abs(pwr_db['clats'] - clat))
                    pwr_lon_ind = 0 # Central longitude only for this version
                    pwr = pwr_db['pwr'][pwr_freq_ind, pwr_lat_ind, pwr_lon_ind]*pow(I0, 2.)
                    # print pwr
                    points_2d = np.hstack([np.vstack([ray_data[k4]['pos'][[0,2],t_ind:t_ind+2], np.zeros([1,2])]),
                                           np.vstack([ray_data[k5]['pos'][[0,2],t_ind:t_ind+2], np.zeros([1,2])]),
                                           np.vstack([ray_data[k6]['pos'][[0,2],t_ind:t_ind+2], np.ones([1,2])*nf]),
                                           np.vstack([ray_data[k7]['pos'][[0,2],t_ind:t_ind+2], np.ones([1,2])*nf])])
                    tri = Delaunay(points_2d.T, qhull_options='QJ')

                    # Mask off the grid points which are definitely outside the volume:
                    minx = min(points_2d[0,:])
                    maxx = max(points_2d[0,:])
                    minz = min(points_2d[1,:])
                    maxz = max(points_2d[1,:])
                    
                    ix = np.where((xx >= minx) & (xx <= maxx))[0]
                    iz = np.where((zz >= minz) & (zz <= maxz))[0]
                    ief= np.arange(0, nf)
                    # print minx, maxx, minz, maxz, len(ix), len(iz)
                    px, pz, pf = np.meshgrid(ix, iz, ief, indexing='ij')  # in 3d, ij gives xyz, xy gives yxz. dumb.
                    newpoints = np.vstack([xx[px.ravel()], zz[pz.ravel()], ff[pf.ravel()]])
                    
                    mask = (tri.find_simplex(newpoints.T) >= 0)*1.0
                    mask = mask.reshape([len(ix), len(iz), len(ief)])
                    total_cells = np.sum(mask)
                    # print total_cells
                    if (total_cells > 0):
                        dc[px,pz,pf] += damping_avg*pwr*mask/voxel_vol
                        # data_cur += np.sum(dc, axis=-1)
                # tmax
            # lat pair
            # data_total[:,:,freq_ind, t_ind] += np.sum(dc, axis=-1)
            data_total[:,:, t_ind] += np.sum(dc, axis=-1)
        # freq pair
    # time step

    outs = dict()

    outs['data'] = data_total
    outs['params'] = dict()
    outs['params']['flash_lat']=flash_lat
    outs['params']['mlt']=mlt
    outs['params']['max_dist']=max_dist
    outs['params']['tmax'] = tmax
    outs['params']['dt']= dt
    outs['params']['I0']=I0
    outs['params']['f_low']=f_low
    outs['params']['f_hi']=f_hi
    outs['params']['step_size'] = step_size
    outs['params']['n_sub_freqs'] = n_sub_freqs
    outs['params']['Llims'] = Llims
    outs['params']['NL'] = NL

    return outs            

# if __name__ == "__main__":

#     logging.basicConfig(level=logging.INFO,
#                         format='[%(levelname)s] %(message)s')  


#     data = interp_ray_power(
#                             ray_dir='/shared/users/asousa/WIPP/rays/2d/nightside/gcpm_kp0',
#                             power_dir = '/shared/users/asousa/WIPP/WIPP_stencils/outputs/input_energies/',
#                             tmax = 10,
#                             flash_lat=40,
#                             mlt=0,
#                             dt=0.05,
#                             f_low=1000,
#                             f_hi=5000,
#                             max_dist=120,
#                             n_sub_freqs=50,
#                             Llims=[0,5],
#                             NL = 100,
#                             dlon = 1,
#                             I0 = 10000
#                             )

#     # # np.save("data_dump.npy",data)
#     # with open('data_dump.pkl','w') as f:
#     #     pickle.dump(data,f)
#     # with gzip.open('data_dump.pklz','wb') as f:
#     #     pickle.dump(data,f)


#   data['data'] has dimensions ([n_fieldlines, n_freq_pairs, n_longitudes, n_times])

    # # Sum over frequencies for plotting:
    # data['data'] = np.sum(data['data'], axis=1)

    # fig, ax = plot_LT(data)
    # fig.savefig("test_timeseries.png",ldpi=300)

    # fig, ax = plot_lon(data)
    # fig.savefig("test_longitude.png",ldpi=300)
    # 