from __future__ import division
from mpi4py import MPI
import commands
from partition import partition 

import matplotlib
matplotlib.use('agg')



import numpy as np
import pandas as pd
import pickle
import gzip

from joblib import Parallel, delayed
import multiprocessing
from index_helpers import load_Kp
import matplotlib.pyplot as plt
import os
import itertools
import random
import os
import time
import datetime as datetime
import types
import scipy.io
import matplotlib.gridspec as gridspec



from scipy import stats
import xflib
import logging
import math
from mpl_toolkits.mplot3d import Axes3D


from mpl_toolkits.axes_grid1 import make_axes_locatable
xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')
from scipy.integrate import nquad

from GLD_file_tools import GLD_file_tools

import aacgmv2

# Input settings:

pwr_db_path = '/shared/users/asousa/WIPP/lightning_power_study/outputs/pwr_db_20deg_spread.pklz'
# output_path = '/shared/users/asousa/WIPP/lightning_power_study/outputs/GLDstats_v6'
output_path = '/shared/users/asousa/WIPP/lightning_power_study/outputs/GLDstats_v9_CGM'  # This one to include Io histogram - 6.29.17
fig_path = os.path.join(output_path, 'figures')
dump_path = os.path.join(output_path,'data')


lookback_time = datetime.timedelta(hours=3)

# The range of times we'll do:
start_day = datetime.datetime(2014,8,9,0,0,0)
stop_day = datetime.datetime(2017,6,1,0,0,0)

CGM = True # Convert from geo to CGM?


# ------------ Start MPI -------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")
nProcs = 1.0*comm.Get_size()


if rank == 0:
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    if not os.path.exists(fig_path):
        os.mkdir(fig_path)
    if not os.path.exists(dump_path):
        os.mkdir(dump_path)
        

# --------------------- Definitions -------------------------------

# Convert the Matlab coastline datafile to geomagnetic coordinates:
def get_coast_mag(itime):
    xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')
    coastlines = scipy.io.loadmat('/shared/users/asousa/WIPP/lightning_power_study/python/coastlines.mat')

    coast_lat_mag = np.zeros(len(coastlines['lat']))
    coast_lon_mag = np.zeros(len(coastlines['long']))

    for ind, (lat, lon) in enumerate(zip(coastlines['lat'], coastlines['long'])):
        if np.isnan(lat) or np.isnan(lon):
            coast_lat_mag[ind] = np.nan
            coast_lon_mag[ind] = np.nan
        else:
            tmpcoords = [1, lat[0], lon[0]]
            tmp_mag = xf.rllgeo2rllmag(tmpcoords, itime)
            coast_lat_mag[ind] = tmp_mag[1]
            coast_lon_mag[ind] = tmp_mag[2]

    # Loop around for -180 + 180 ranges
    coast_lat_mag = np.concatenate([coast_lat_mag, coast_lat_mag[coast_lon_mag > 180]])
    coast_lon_mag = np.concatenate([coast_lon_mag, (coast_lon_mag[coast_lon_mag > 180] - 360)])

    # Toss in some NaNs to break up the continents
    for ind in range(len(coast_lat_mag) -1):
        if ((np.abs(coast_lat_mag[ind+1] - coast_lat_mag[ind]) > 5) or
           (np.abs(coast_lon_mag[ind+1] - coast_lon_mag[ind]) > 5)):
            coast_lat_mag[ind] = np.nan
            coast_lon_mag[ind] = np.nan

    return coast_lat_mag, coast_lon_mag

def data_grid_at(in_time):
    # P_A = 5e3
    # P_B = 1e5
    # tpeak  = np.log(P_A/P_B)/(P_A - P_B)
    # Ipeak = np.exp(-P_A*tpeak) - np.exp(-P_B*tpeak)
    # Ipeak2Io = 1.0/Ipeak
    Ipeak2Io = 1.2324    # Conversion between peak current and Io
                         # (i.e., normalize the exponential terms)

#     print np.shape(times_to_do)
    print "loading flashes at ", in_time
    data_grid = []
#     for in_time in times_to_do:

    Kpm = Kpmax[Kpmtimes == in_time]
    Kp_cur = Kp[np.where(Ktimes == in_time)[0]]

    flashes, flash_times = gld.load_flashes(in_time, lookback_time)
    print np.shape(flashes)
    if flashes is not None:
        if CGM:
            cgmlat, cgmlon = aacgmv2.convert(flashes[:,7], flashes[:,8], 5.0*np.ones_like(flashes[:,7]))
            happy_inds = ~np.isnan(cgmlat)
            cgmlat = cgmlat[happy_inds]
            cgmlon = cgmlon[happy_inds]
            cgmlon[cgmlon < 0] += 360.
            mlts = [xf.lon2MLT(ft, cl) for (ft, cl) in zip(flash_times[happy_inds], cgmlon)] # Eh, good enough. (This really should be done on unconverted dipole instead of CGM)

            I = flashes[happy_inds,9]*Ipeak2Io

            data_grid = np.vstack([cgmlat, cgmlon, mlts, I, Kpm*np.ones_like(I), Kp_cur*np.ones_like(I)]).T
        else:
            for flash, flashtime in zip(flashes, flash_times):

                glat = flash[7]
                glon = flash[8]
                # I    = flash[9]*Ipeak2Io  # Added after stats_v6
                I    = flash[9]*Ipeak2Io  # Added after stats_v6
                # Get location in geomagnetic coordinates
                mloc = xf.rllgeo2rllmag([1.0, glat, glon], flashtime)

                # Get MLT:
                mlt = xf.lon2MLT(flashtime, mloc[2])

                data_grid.append([mloc[1], mloc[2], mlt, I, Kpm, Kp_cur])
            data_grid = np.array(data_grid)

        return data_grid
    else:
        return None


# Output power space:
def analyze_flashes(data_grid, in_time):
    Ipeak2Io = 1.2324    # Conversion between peak current and Io

    outdata = dict()
    dg2 = np.array(data_grid)

    # Quantize data into lat, lon, and MLT bins:
    dg2[:,0] = np.digitize(dg2[:,0],gridlats)
    dg2[:,1] = np.digitize(np.mod(dg2[:,1], 360.0) - 180.0, gridlons)
    dg2[:,2] = (dg2[:,2] > 6) & (dg2[:,2] <= 18)   # Is day?

    hist_bins = np.arange(0,24.5, 0.5)
    flash_map = np.zeros([len(gridlats), len(gridlons)])
    cur_map   = np.zeros([len(gridlats), len(gridlons)])
    pwr_map   = np.zeros([len(gridlats), len(gridlons)])
    mlt_hist, _  = np.histogram(data_grid[:,2], hist_bins)

    # Io_bins = np.linspace(0,1000,101)
    Io_bins = np.unique(np.round(pow(10,np.linspace(0,3,144))))*Ipeak2Io # 101 log-spaced dividers
    Io_hist, _ = np.histogram(np.abs(data_grid[:,3]), Io_bins)

    # Bin total current by lat and lon
    day_bins = np.zeros([len(gridlats), len(gridlons)])
    nite_bins= np.zeros([len(gridlats), len(gridlons)])

    
    for row in dg2:
        # Power stencils are proportional to current squared;
        # GLD currents are in kA
        if row[2]:
            day_bins[int(row[0]), np.mod(int(row[1]), 360)] += pow(row[3]*1e3, 2.0)
        else:
            nite_bins[int(row[0]), np.mod(int(row[1]), 360)] += pow(row[3]*1e3, 2.0)
        
        # Flash count histogram
        flash_map[int(row[0]), np.mod(int(row[1]), 360)] += 1
        
        # 2d current histogram
        cur_map[int(row[0]), np.mod(int(row[1]), 360)] += pow(row[3]*1e3, 2.0)


# --------- Commenting out the power stencils for this version
    # day_todo = np.where(day_bins > 0)
    # nite_todo = np.where(nite_bins > 0)

    # for isday in [False, True]:
    #     if isday:
    #         todo = np.where(day_bins > 0)
    #     else:
    #         todo = np.where(nite_bins > 0)
    #     for latind, lonind in zip(todo[0], todo[1]):
    #         if (np.abs(gridlats[latind]) >= stencil_lats[0]) & (np.abs(gridlats[latind]) <= stencil_lats[-1]):

    #             if isday:
    #                 key = (np.round(gridlats[latind]), 12)
    #             else:
    #                 key = (np.round(gridlats[latind]), 0)
    #             if key in inp_pwr_dict:
    #                 stencil = inp_pwr_dict[key]
    #                 if isday:
    #                     pwr = day_bins[latind, lonind]
    #                 else:
    #                     pwr = nite_bins[latind, lonind]
    #                 latleft = int(latind - cell_lat_offset)
    #                 latright = int(latind + cell_lat_offset-1)
    #                 lonleft = int(lonind - cell_lon_offset)
    #                 lonright =int(lonind + cell_lon_offset-1)

    #                 if lonleft < 0:
    #                     # Wrap around left:
    #                     pwr_map[latleft:latright, 0:lonright] += \
    #                             stencil[:, np.abs(lonleft):]*pwr
    #                     pwr_map[latleft:latright, (len(gridlons) - np.abs(lonleft)):] += \
    #                             stencil[:,0:np.abs(lonleft)]*pwr
    #                 elif lonright >= len(gridlons):
    #                     # wrap around right:
    #                     pwr_map[latleft:latright, lonleft:len(gridlons)] += \
    #                         stencil[:,0:len(gridlons) - lonleft]*pwr
    #                     pwr_map[latleft:latright, 0:np.abs(lonright) - len(gridlons)] += \
    #                         stencil[:,len(gridlons) - lonleft:]*pwr
    #                 else:
    #                     pwr_map[latleft:latright, lonleft:lonright] += stencil*pwr
    # -----

    # Roll the output data arrays to get [-180, 180] instead of [0, 360]
    # outdata['pwr_map']   = np.roll(pwr_map,  len(gridlons)/2, axis=1)
    outdata['flash_map'] = np.roll(flash_map,int(len(gridlons)/2), axis=1)
    outdata['cur_map']   = np.roll(cur_map , int(len(gridlons)/2), axis=1)
    outdata['mlt_hist']  = mlt_hist
    outdata['Io_hist']   = Io_hist
    outdata['in_time']   = in_time
    return outdata

def plot_pwr_data(outdata):
    # --------------- Latex Plot Beautification --------------------------
    fig_width = 12 
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
    
    pwr_map = outdata['pwr_map']
    flash_map = outdata['flash_map']
    cur_map = outdata['cur_map']
    mlt_hist = outdata['mlt_hist']
    in_time = outdata['in_time']

    clims=[0, 4] # log space

    # Plot flash location with integrated current vs lat and lon:
    fig = plt.figure()
    gs = gridspec.GridSpec(2,3,width_ratios=[0.5,10,0.25],
                       height_ratios=[10,1])
    # gs.update(left=0.05, right=0.48, wspace=0.05)
    ax1 = plt.subplot(gs[0:-1,1:-1])  # main figure
    ax0 = plt.subplot(gs[0:-1,0])   #
    ax2 = plt.subplot(gs[-1,1:-1])
    cbar_ax = plt.subplot(gs[:, -1])
    pwr_bylat = np.sum(pwr_map, axis=1)/lookback_time.total_seconds()
    pwr_bylon = np.sum(pwr_map, axis=0)/lookback_time.total_seconds()
    
    logpwr = np.log10(pwr_map/3600./3.0)  # Energy per second, in log space.
    logpwr[np.isinf(logpwr)] = -10


    p = ax1.pcolorfast(gridlons, gridlats, logpwr, cmap = plt.get_cmap('viridis'))
    p.set_clim(clims)

    cb = plt.colorbar(p, cax=cbar_ax)

    cb.set_label('Average Energy Flux [J/sec]')
    cticks = np.arange(clims[0],clims[1] + 1)
    cb.set_ticks(cticks)
    cticklabels = ['$10^{%d}$'%k for k in cticks]
    cb.set_ticklabels(cticklabels)



    coast_lat_mag, coast_lon_mag = get_coast_mag(in_time)
    ax1.plot(coast_lon_mag, coast_lat_mag, 'w')
    ax1.set_xlim([-180, 180])
    ax1.set_ylim([-90, 90])
    p.set_clim([0,4])
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax0.plot(pwr_bylat,gridlats)
    ax0.set_ylim([-90,90])
    ax0.set_xticks([])
    ax0.set_ylabel('Latitude (magnetic)')
    ax2.plot(gridlons, pwr_bylon)
    ax2.set_xlim([-180,180])
    ax2.set_yticks([])
    ax2.set_xlabel('Longitude (magnetic)')

    ax1.set_title( in_time.isoformat() )
    fig.tight_layout()

    return fig, ax0, ax1, ax2

# ----------------- MAIN RUN BLOCK --------------------------
# Load input power database:
with gzip.open(pwr_db_path,'rb') as f:
    inp_pwr_dict = pickle.load(f)

tmp = np.array(sorted([k for k in inp_pwr_dict.keys() if not isinstance(k, basestring)]))

# Flip the stencils for southern hemisphere
for k in tmp:
    stencil = inp_pwr_dict[tuple(k)]
    newkey = (-1*k[0], k[1])
    inp_pwr_dict[newkey] = np.flipud(stencil)

tmp = np.array(sorted([k for k in inp_pwr_dict.keys() if not isinstance(k, basestring)]))
stencil_lats = np.unique(tmp[:,0]); stencil_MLTs = np.unique(tmp[:,1])

cellsize = inp_pwr_dict['cellsize']
cell_lat_offset = np.round(inp_pwr_dict['lat_spread']/inp_pwr_dict['cellsize'])
cell_lon_offset = np.round(inp_pwr_dict['lon_spread']/inp_pwr_dict['cellsize'])

gridlats = np.arange(-90, 90, cellsize)
gridlons = np.arange(-180, 180, cellsize)
# gridmlts = np.linspace(0,24, 9)

# Get Kp data
Ktimes, Kp = load_Kp()
Ktimes = [k + datetime.timedelta(minutes=90) for k in Ktimes]  # 3-hour bins; the original script labeled them in the middle of the bin
Ktimes = np.array(Ktimes)
Kp = np.array(Kp)

# Get Kpmax -- max value of Kp over the last 24 hours (8 bins):
Kpmax = np.max([Kp[0:-8],Kp[1:-7],Kp[2:-6], Kp[3:-5], Kp[4:-4],Kp[5:-3],Kp[6:-2], Kp[7:-1], Kp[8:]],axis=0)
Kpmtimes = Ktimes[8:]


# lookback_time = np.unique(np.diff(Kpmtimes))[0] # This should be 3 hours

# num_cores = multiprocessing.cpu_count()


# Lightning data getter
GLD_path = '/home/asousa/GLD_mount/'
gld = GLD_file_tools(GLD_path, prefix='GLD')



# Do one file per day, since we're just parallelizing over cores.
# 8 kp per day, 8 cores, bam.
# def job(intime):
#     datagrid = data_grid_at(intime)
#     if datagrid is not None:
#         return analyze_flashes(datagrid, intime)
#     else:
#         return None



# days_to_do = [start_day + datetime.timedelta(days=x) for x in range((stop_day - start_day).days + 1)]
# day_pairs = zip(days_to_do[0:-1], days_to_do[1:])

# for day1, day2 in day_pairs:
#     times_to_do = Ktimes[(Ktimes > day1) & (Ktimes <= day2)]
#     # print times_to_do
#     datae = Parallel(n_jobs=num_cores)(delayed(job)(t) for t in times_to_do)

#     for datum in datae:
#         if datum is not None:
#             fig, ax0, ax1, ax2 = plot_pwr_data(datum)
#             figname = os.path.join(fig_path, datum['in_time'].strftime('%d_%m_%Y_%H_%M') +'.png')
#             fig.savefig(figname, ldpi=300)
#             plt.close('all')

#             filename = os.path.join(dump_path, datum['in_time'].strftime('%d_%m_%Y_%H_%M') + '.pklz')
#             with gzip.open(filename, 'wb') as f:
#                 pickle.dump(datum,f)

# days_to_do = [start_day + datetime.timedelta(days=x) for x in range((stop_day - start_day).days + 1)]
# day_pairs = zip(days_to_do[0:-1], days_to_do[1:])


if rank == 0:
    tasklist = Ktimes[(Ktimes > start_day) & (Ktimes <= stop_day)]
    chunks = partition(tasklist, nProcs)
else:
    tasklist = None
    chunks = None

tasklist = comm.bcast(tasklist, root=0)
chunks   = comm.bcast(chunks, root=0)

nTasks  = 1.0*len(tasklist)
nSteps = np.ceil(nTasks/nProcs).astype(int)

if (rank < len(chunks)):
    print "Subprocess %s on %s: doing %d tasks"%(rank, host, len(chunks[rank]))


    for intime in chunks[rank]:
        datagrid = data_grid_at(intime)
        if datagrid is not None:

            filename = os.path.join(dump_path, intime.strftime('%m_%d_%Y_%H_%M') + '.pklz')
            datum = analyze_flashes(datagrid, intime)

            if datum is not None:
                # fig, ax0, ax1, ax2 = plot_pwr_data(datum)
                # figname = os.path.join(fig_path, datum['in_time'].strftime('%m_%d_%Y_%H_%M') +'.png')
                # fig.savefig(figname, ldpi=300)
                # plt.close(fig)

                filename = os.path.join(dump_path, datum['in_time'].strftime('%m_%d_%Y_%H_%M') + '.pklz')
                with gzip.open(filename, 'wb') as f:
                    pickle.dump(datum,f)




#     datae = Parallel(n_jobs=num_cores)(delayed(job)(t) for t in times_to_do)

#     for datum in datae:
#         if datum is not None:
#             fig, ax0, ax1, ax2 = plot_pwr_data(datum)
#             figname = os.path.join(fig_path, datum['in_time'].strftime('%d_%m_%Y_%H_%M') +'.png')
#             fig.savefig(figname, ldpi=300)
#             plt.close('all')

#             filename = os.path.join(dump_path, datum['in_time'].strftime('%d_%m_%Y_%H_%M') + '.pklz')
#             with gzip.open(filename, 'wb') as f:
#                 pickle.dump(datum,f)



