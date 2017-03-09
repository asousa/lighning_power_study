import matplotlib
matplotlib.use('agg')
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

def interp_ray_power(ray_dir='/shared/users/asousa/WIPP/rays/2d/nightside/gcpm_kp0/',
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
                    d_lon = 1,
                    num_lons = 10,
                    Llims = [1.2, 5],
                    L_step = 0.05,
                    frame_directory=None):

    # Constants
    Hz2Rad = 2.*np.pi
    D2R = np.pi/180.
    H_IONO_BOTTOM = 1e5
    H_IONO_TOP = 1e6
    R_E = 6371e3


   
    Lshells = np.arange(Llims[0], Llims[1], L_step)
    L_MARGIN = L_step/2.0



    # Coordinate transform tools
    xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')

    t = np.arange(0,tmax, dt)
    itime = datetime.datetime(2010,1,1,0,0,0)

    # Select flash longitude corresponding to the input MLT:
    flash_lon = MLT2lon(itime,mlt,xf)

    # if ((mlt >= 0) and (mlt <= 6)) or ((mlt >= 18) and (mlt <= 24)):
    #     sun = xf.gse2sm([-1,0,0], itime)
    #     sun_geomag_midnight = (xf.sm2rllmag(sun, itime))
    #     flash_lon = sun_geomag_midnight[2]
    # else:
    #     sun = xf.gse2sm([1,0,0], itime)
    #     sun_geomag_noon = (xf.sm2rllmag(sun, itime))
    #     flash_lon = sun_geomag_noon[2]

    print "flash MLT:",mlt, "Flash longitude:",flash_lon
    flash_lons = np.arange(flash_lon, flash_lon + (num_lons +1)*d_lon, d_lon)


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

    newlons = np.array([flash_lon - dlat/2., flash_lon + dlat/2.])
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
        logging.info("Loading freq: %d"%freq)
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
    logging.info("Rotating to new longitudes")
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

    # print raylats, raylons

    inp_pwrs = dict()

    lat_pairs  = zip(raylats[0:-1],raylats[1:])
    lon_pairs  = zip(flash_lons[0:-1],flash_lons[1:])
    # print lon_pairs

    flash_pos_mag = [1, flash_lat, flash_lon]
    flash_pos_sm = xf.rllmag2sm(flash_pos_mag, itime)


    opts = dict()
    opts['epsabs']= 1.5e-8
    opts['epsrel']= 1.5e-8
    opts['limit']= 10

    def integrand(inlat, inlon, inw, itime, I0, flash_pos_sm_in, itime_in):
        mlt = lon2MLT(itime, inlon, xf);
        # print "lon:", inlon, "MLT:",mlt
        tmp_coords = [1, inlat, inlon];
        x_sm = xf.rllmag2sm(tmp_coords, itime_in);

        pwr = input_power_scaling(flash_pos_sm_in, x_sm, inlat, inw, I0, mlt, xf);
        return pwr*(R_E + H_IONO_TOP)*D2R*(R_E + H_IONO_TOP)*np.cos(D2R*inlat)*D2R

    # pwr_vec = np.zeros(len(lat_pairs))
    inp_pwrs = dict()
    for f1, f2 in freq_pairs:
        logging.info("\t%d, %d"%(f1, f2))
        # n_freqs = np.ceil(np.abs(f2 - f1)/min_fstep)
        # f_weights = (np.arange(0,1,1.0/n_sub_freqs) + (1.0/(2.*n_sub_freqs)))
        
        # for f_weight in f_weights:
        # f_center = f_weight*f1 + (1.0-f_weight)*f2
        f_center = (f1 + f2)/2.

        for ind, (lat1, lat2) in enumerate(lat_pairs):
            pwr_vec = np.zeros(len(lon_pairs))

            clat = (lat1 + lat2)/2.
            w1 = Hz2Rad*f1
            w2 = Hz2Rad*f2
            w   = Hz2Rad*(f1 + f2)/2.
            dw = np.abs(f1 - f2)*Hz2Rad
            for lon_ind, lon_pair in enumerate(lon_pairs):
                ranges = [[lat1, lat2], lon_pair]
                # Integrate power in latitude and longitude
                integ = nquad(integrand, ranges, args=[w, itime, I0, flash_pos_sm, itime], opts=opts, full_output=False)
                pwr = integ[0]
                pwr_vec[lon_ind] = pwr*dw

            key = (f_center, clat)
            inp_pwrs[key] = pwr_vec

    logging.info('Total input energy: %0.1f J'%(np.sum(inp_pwrs.values())))

    with open('input_energy_%d.pkl'%flash_lat,'w') as file:
        pickle.dump(inp_pwrs, file)




# ------------ Set up L-shell output grid --------------------
# Similar deal, but now let's look along field lines instead of a Cartesian grid

    # n_lsteps = 100
    R2D = 180./np.pi
    D2R = np.pi/180.

# ------------------ Set up field lines ----------------------------
    fieldlines = []
    for L in Lshells:
        fieldline = dict()
        maxlat = np.floor(np.arccos(np.sqrt((R_E + H_IONO_TOP)/R_E/L))*R2D)
        n_lsteps = int(np.round(2.0*maxlat/dlat))
        lat_divisions = np.linspace(maxlat, -1.0*maxlat, n_lsteps+1)
        lat_centers   = lat_divisions[0:-1] - dlat/2.
    
        fieldline['lat'] = lat_centers
        
        # Radius of tube around field line:
        clam = np.sin(lat_centers*D2R)
        slam = np.cos(lat_centers*D2R)
        clam2 = pow(clam,2.)
        slam2 = pow(slam,2.)
        rootTerm = np.sqrt(1.0*3.0*slam2)
       
        clam = np.cos(lat_centers*D2R);
        slam = np.sin(lat_centers*D2R);
        clam2 = pow(clam,2);
        slam2 = pow(slam,2);
        rootTerm = np.sqrt(1+3*slam2);

        radii = clam2*clam / rootTerm * L_MARGIN
        R_centers = L*clam2
        
        fieldline['R'] = R_centers
        fieldline['xradius']= radii
        
        # Approximate each segment as a cylinder:
        seg_length = R_centers*dlat*D2R
        seg_vol = np.pi*pow(radii,2.)*seg_length*pow(R_E,3.)  # cubic meters

        fieldline['vol'] = seg_vol
        fieldline['total_vol'] = np.sum(seg_vol)

        fieldline['x'] = R_centers*clam
        fieldline['y'] = R_centers*slam
        
        fieldline['x_unit_vect'] = (3*clam2 - 2) / rootTerm ;
        fieldline['y_unit_vect'] = (3*slam*clam) / rootTerm ;

        coords_rllmag = np.vstack([R_centers, lat_centers, np.ones(len(lat_centers))*flash_lon])
        coords_sm = []
        for row in coords_rllmag.T:
            coords_sm.append(xf.rllmag2sm(row, itime))
        fieldline['pos'] = np.array(coords_sm)

        fieldlines.append(fieldline)

#----------- Step through and fill in the voxels (the main event) ---------------------
    logging.info("Starting interpolation")

    # output space
    nfl = len(fieldlines)
    nlons = len(flash_lons) - 1
    nt = len(t)
    n_freq_pairs = len(freq_pairs)
    data_total = np.zeros([nfl, n_freq_pairs, nlons, nt])

    lon1 = raylons[0]
    lon2 = raylons[1]

    for t_ind in np.arange(nt - 1):
        # Per frequency
        data_cur = np.zeros(nfl)
        logging.info("t = %d"%t_ind)
        for freq_ind, (f1, f2) in enumerate(freq_pairs):
            # Loop over adjacent sets:
            ff = np.arange(0, n_sub_freqs, 1)
            nf = len(ff)

            for lat1, lat2 in lat_pairs:
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
                    
                    pwr = inp_pwrs[(f_center, clat)]  # Avg pwrs per frequency bin
                    
                    points_2d = np.hstack([np.vstack([ray_data[k4]['pos'][[0,2],t_ind:t_ind+2], np.zeros([1,2])]),
                                           np.vstack([ray_data[k5]['pos'][[0,2],t_ind:t_ind+2], np.zeros([1,2])]),
                                           np.vstack([ray_data[k6]['pos'][[0,2],t_ind:t_ind+2], np.ones([1,2])*nf]),
                                           np.vstack([ray_data[k7]['pos'][[0,2],t_ind:t_ind+2], np.ones([1,2])*nf])])
                    tri = Delaunay(points_2d.T, qhull_options='QJ')
                    
                    for fl_ind, fl in enumerate(fieldlines):
                        ix = np.arange(0,len(fl['pos']))
                        ief= np.arange(0, nf)    
                        px, pf = np.meshgrid(ix, ief, indexing='ij')  # in 3d, ij gives xyz, xy gives yxz. dumb.
                        newpoints = np.hstack([fl['pos'][px.ravel(),:][:,[0,2]], np.atleast_2d(ff[pf.ravel()]).T])

                        mask = (tri.find_simplex(newpoints) >= 0)*1.0
                        mask = mask.reshape([len(ix), len(ief)])
                        
                        total_cells = np.sum(mask)
                        if (total_cells > 0):
                            data_total[fl_ind, freq_ind, :, t_ind] += (damping_avg*pwr/voxel_vol)* \
                                            np.sum(mask*fl['vol'][:,np.newaxis])/fl['total_vol']
                        
        # print "Energy: ", np.sum(data_total[:, :, t_ind], axis=0)

    logging.info("finished with interpolation")

    out_data = dict()
    out_data['data'] = data_total
    out_data['time'] = t
    out_data['Lshells'] = Lshells
    out_data['lons'] = flash_lons
    out_data['flash_lat'] = flash_lat
    out_data['I0'] = I0
    out_data['fmin'] = f_low
    out_data['fmax'] = f_hi
    out_data['freq_pairs'] = freq_pairs


    return out_data


# Plot energy density per L-shell vs time:
def plot_LT(data_dict, tlims=None, Llims=None, clims=None):
    data = data_dict['data']
    Lvec = data_dict['Lshells']
    tvec = data_dict['time']
#     nL, nLons, nt = np.shape(data)
    
#     Lvec = np.arange(Llims[0], Llims[1], dl)
#     tvec = np.arange(0, nt, 1)*dt
    
    fig, ax = plt.subplots(1,1)
    logdata = np.log10(data[:,0,:])
    logdata[np.isinf(logdata)] = -100

    if clims is None:
        maxval = np.ceil(np.max(logdata))
        clims = [maxval - 5, maxval]
    
    print np.shape(data)
    p0 = ax.pcolorfast(tvec, Lvec, logdata, vmin=clims[0], vmax=clims[1], cmap=plt.get_cmap('jet'))

    fig.subplots_adjust(right=0.82)
    cax = fig.add_axes([0.84,0.12, 0.02, 0.75])

    cb = plt.colorbar(p0, cax=cax)
    cb.set_label('avg wave power density (J/m$^3$)')
    cticks = np.arange(clims[0],clims[1] + 1)
    cb.set_ticks(cticks)
    cticklabels = ['$10^{%d}$'%k for k in cticks]
    cb.set_ticklabels(cticklabels)



    ax.set_ylabel('L shell')
    ax.set_xlabel('Time (sec)')
    if tlims is not None:
        ax.set_xlim(tlims)

    return fig, ax
    
# plot_LT(data_total,[0, 20],dt, [1.2, 5], L_step)

# 2d plot, longitude vs L-shell:
def plot_lon(data, lonlims=None, Llims=None, clims=None):
    # --------------- Latex Plot Beautification --------------------------
    fig_width = 8 
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
    
    data_total = data['data']
    flash_lons = data['lons']
    Lshells    = data['Lshells']
    flash_lon  = flash_lons[0]
    time = data['time']
    dt = time[1]-time[0]
    tmax= time[-1] + dt
    d_lon = flash_lons[1] - flash_lons[0]
    num_lons = len(flash_lons)

    lon_max = flash_lons[-1] - flash_lons[0]
        
    logdata = np.log10(np.sum(data_total, axis=-1)*dt/tmax)
    logdata[np.isinf(logdata)] = -1000
    logdata = np.vstack([np.flipud(logdata[:,1:].T), logdata.T]).T

    maxval = np.ceil(np.max(logdata))
    clims = [maxval - 3, maxval]

    fig, ax = plt.subplots(1,1)
    xaxis = np.arange(-lon_max, lon_max + d_lon, d_lon)
    p0 = ax.pcolorfast(xaxis, Lshells, logdata, vmin=clims[0], vmax=clims[1], cmap=plt.get_cmap('jet'))

    fig.subplots_adjust(right=0.82)
    cax = fig.add_axes([0.84,0.12, 0.02, 0.75])

    cb = plt.colorbar(p0, cax=cax)
    cb.set_label('avg wave power density (J/sec/m$^3$)')
    cticks = np.arange(clims[0],clims[1] + 1)
    cb.set_ticks(cticks)
    cticklabels = ['$10^{%d}$'%k for k in cticks]
    cb.set_ticklabels(cticklabels)

    ax.set_xlabel('Longitude (deg)')
    ax.set_ylabel('L shell')

    return fig, ax






if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO,
                        format='[%(levelname)s] %(message)s')  

    Llims = [1.2, 6]
    L_step = 0.05

    d_lon = 0.5
    num_lons = 100

    data = interp_ray_power(ray_dir='/shared/users/asousa/WIPP/rays/2d/nightside/gcpm_kp0',
                                tmax = 10,
                                flash_lat=40,
                                mlt=0,
                                dt=0.1,
                                f_low=10000,
                                f_hi=30000,
                                max_dist=300,
                                n_sub_freqs=50,
                                Llims=Llims,
                                L_step=L_step,
                                d_lon = d_lon,
                                num_lons = num_lons
                                )

    # np.save("data_dump.npy",data)
    with open('data_dump.pkl','w') as f:
        pickle.dump(data,f)
    with gzip.open('data_dump.pklz','wb') as f:
        pickle.dump(data,f)


#   data['data'] has dimensions ([n_fieldlines, n_freq_pairs, n_longitudes, n_times])

    # Sum over frequencies for plotting:
    data['data'] = np.sum(data['data'], axis=1)

    fig, ax = plot_LT(data)
    fig.savefig("test_timeseries.png",ldpi=300)

    fig, ax = plot_lon(data)
    fig.savefig("test_longitude.png",ldpi=300)
    




