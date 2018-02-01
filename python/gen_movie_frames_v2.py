import numpy as np
import pandas as pd
import pickle
#from build_database import flux_obj
from scipy import interpolate
import matplotlib.pyplot as plt
import os
import itertools
import random
import os
import time
import datetime
import gzip
import matplotlib.gridspec as gridspec

import xflib as xflib

R2D = 180./np.pi
D2R = np.pi/180.


from ray_interp_2d_grid_v2 import interp_ray_power

def plot_frame(data, t_ind, clims=None):
    # --------------- Latex Plot Beautification --------------------------
    fig_width = 7 
    fig_height = 5
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
    Kp = 0
    Lpp,_,_ = bulge(data['params']['mlt'], Kp)
    print Lpp
    ninterp = 800
    cm = plt.get_cmap('jet')
    
    Llims = data['params']['Llims']
    NL = data['params']['NL']
    
    xx = np.linspace(Llims[0], Llims[1], NL)
    zz = np.linspace(-Llims[1]/2., Llims[1]/2., NL)
    
    newx = np.linspace(Llims[0], Llims[1], ninterp)
    newz = np.linspace(-Llims[1]/2., Llims[1]/2., ninterp)

    newx = np.linspace(0, 6, ninterp)
    newz = np.linspace(-3, 3, ninterp)


#     Originally: nx, nz, nf, nt. 
#     Select time frame and um over freq axis for this plot
#     d = np.sum(data['data'][:,:,:,t], axis=-1)

    # nx, nx, nt. Select current timestep.
    d = data['data'][:,:,t]
#     logdata = np.log10(d)
#     logdata[np.isinf(logdata)] = -100


    px, pz = np.meshgrid(xx, zz)
    
    
    interp = interpolate.RegularGridInterpolator([xx, zz], d, bounds_error=False, fill_value=0)
    px, py = np.meshgrid(newx, newz)
    pts = zip(px.ravel(), py.ravel())
    logdata = np.log10(interp(pts)).reshape(len(newx), len(newz))
    logdata[np.isinf(logdata)] = -100

    # Load coastlines (for plotting)
    with gzip.open('../../Thesis figures/python_local/mag_coastlines.gzip','rb') as file:
        coast = pickle.load(file)
    
    maxlog = np.max([logdata])
    
    # Show about 5 orders of magnitude
    if clims is None:
        clims = [maxlog - 5, maxlog]

#     fig, ax = plt.subplots(1,1)
    fig = plt.figure()
    gs = gridspec.GridSpec(1,2, width_ratios=[1,0.05])
    gs.update(wspace=0.05, hspace=0.1) # set the spacing between axes.
    
    ax = plt.subplot(gs[0])
    cax = plt.subplot(gs[1])

    # Plot the earth
    earth = plt.Circle((0,0),1,color='0.5',alpha=1, linewidth=2,zorder=100)
    iono  = plt.Circle((0,0),(R_E + H_IONO_TOP)/R_E, color='w',alpha=0.8, zorder=99)
    ax.add_patch(earth)   
    ax.add_patch(iono)
        
    # Plot coastlines
    coastpoints = np.vstack([coast['lon']/90. + 0.3, coast['lat']/90.])
    coastr = np.linalg.norm(coastpoints, axis=0)
    coastmask = (coastr < 1) | (np.isnan(coastr))
    ax.plot(coastpoints[0,coastmask], coastpoints[1,coastmask],'k', zorder=101, alpha=0.8, linewidth=1)

    # Plot plasmapause fieldline:
    
    flats = np.arange(-80,80,1)*D2R
    fl_rad = Lpp*pow(np.cos(flats),2)
    ax.plot(fl_rad*np.cos(flats), fl_rad*np.sin(flats),'--',color='w',linewidth=2)
    
    # Plot data!
    p0 = ax.pcolorfast(newx, newz, logdata, vmin=clims[0], vmax=clims[1], cmap = cm)
    ax.set_aspect('equal')
    ax.set_xlim([newx[0],newx[-1]])
    ax.set_ylim([newz[0],newz[-1]])

    # Colorbar
    cb = plt.colorbar(p0, cax=cax)
    cticks = np.arange(clims[0],clims[1] + 1)
    cb.set_ticks(cticks)
    cticklabels = ['$10^{%d}$'%k for k in cticks]
    cb.set_ticklabels(cticklabels)
    cb.set_label('Energy density [J/m$^3$]') 
    
    fig.subplots_adjust(right=0.8)

    return fig









flash_lat = int(os.getenv('inlat'))
kp   = int(os.getenv('kp'))
outdir = '/shared/users/asousa/WIPP/lightning_power_study/outputs/movie_frames_jan_2018/inlat_%d'%flash_lat

if not os.path.exists(outdir):
    os.system('mkdir -p %s'%outdir) 

print "doing thing for kp = ", kp
data = interp_ray_power(
        ray_dir='/shared/users/asousa/WIPP/rays/2d/nightside/mode6/kp%d/'%kp,
        power_dir = '/shared/users/asousa/WIPP//Input_Power_Debugging_2018/outputs/input_energies_1lax0.25lox33f_one_sided/',
        tmax = 20,
        flash_lat=flash_lat,
        mlt=0,
        dt=0.04,
        f_low=200,
        f_hi=30000,
        max_dist=1200,
        n_sub_freqs=50,
        Llims=[0,8],
        NL = 400,
        dlon = 0.25,
        I0 = 10000
        )

data['params']['kp'] = kp

print "Saving for kp%d"%kp
# with gzip.open(os.path.join(outdir, 'movie_frames_kp%d.pklz'%kp),'wb') as file:
with open(os.path.join(outdir, 'movie_frames_kp%d.pkl'%kp),'wb') as file:
    pickle.dump(data,file, pickle.HIGHEST_PROTOCOL)






