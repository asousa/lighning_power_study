import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


from mpi4py import MPI
import numpy as np
from scipy import interpolate
import os
from partition import partition
import itertools
import random
import time
import datetime as dt
import commands
import subprocess
import shutil

from index_helpers import load_TS_params
from index_helpers import load_Dst
from index_helpers import load_Kp
from index_helpers import load_ae
from index_helpers import Kp_at
from index_helpers import Ae_at

from raytracer_utils import readdump


import xflib  # Fortran xform-double library (coordinate transforms)
import bisect   # fast searching of sorted lists
from bmodel_dipole import bmodel_dipole
# ------------ Start MPI -------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")
print host
if "cluster" in host:
    host_num = int(host[-3:])
elif "nansen" in host:
    host_num = 0


print "host: %s, host num: %d"%(host, host_num)
nProcs = 1.0*comm.Get_size()


# ------------------Constants --------------------------
R_E = 6371.0
H_IONO=1000.0


# -------------- Simulation params ---------------------
modelnum = 2    # Which model to use (1 = ngo, 2=GCPM, 3=GCPM interp, 4=GCPM rand interp)
use_IGRF = 1    # Magnetic field model (1 for IGRF, 0 for dipole)
use_tsyg = 0    # Use the Tsyganenko magnetic field model corrections

minalt   = (R_E + 150)*1e3 # cutoff threshold in meters

dump_model = True


maxD = 10.0
clims = [-1, 5]

# ---------------- Input parameters --------------------

sim_start = dt.datetime(2001, 01, 1, 0, 0, 0)
sim_stop  = dt.datetime(2002, 01, 1, 0, 0, 0)
sim_step  = dt.timedelta(days=1)


project_root = '/shared/users/asousa/WIPP/lightning_power_study/'
raytracer_root = '/shared/users/asousa/software/raytracer_v1.17/'
ray_bin_dir    = os.path.join(raytracer_root, 'bin')
ray_out_dir = os.path.join(project_root, 'rays','model_dumps_2001')

# GCPM grid to use (plasmasphere model)
if modelnum==1:
    configfile = os.path.join(project_root, 'python','newray.in')
# interpfile = os.path.join(project_root,'raytracer_runscripts','gcpm_models','gcpm_kp40_20010101_0000_MLD01.txt')
if modelnum==3:
    interpfile = os.path.join(project_root,'gcpm_models','demo_models','gcpm_kp4_2001001_L10_80x80x80_noderiv.txt')
if modelnum==4:
    interpfile = os.path.join(project_root, 'gcpm_models','demo_models','gcpm_kp4_2001001_L10_random_5000_20000_0_200000_600000.txt')
    scattered_interp_window_scale = 1.5
    scattered_interp_order = 2
    scattered_interp_exact = 0
    scattered_interp_local_window_scale = 5



timelist = []
if rank == 0:
    curtime = sim_start
    while curtime < sim_stop:
        timelist.append(curtime)
        curtime += sim_step
comm.bcast(timelist, root=0)

planes = ['XY', 'XZ', 'YZ']


# Coordinate transformation library
xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')

# Sync up:
comm.Barrier()

working_path = os.path.join(os.path.expanduser("~"),"rayTmp")

# (if using full GCPM model, you need all the stupid data files in your working directory)
os.chdir(working_path)

# Dump plasmasphere model
if (dump_model):

    if rank==0:
        if not os.path.exists(ray_out_dir):
            os.mkdir(ray_out_dir)


        # tasklist = [(d,p) for d in timelist for p in planes]
        tasklist = timelist
        chunks = partition(tasklist, nProcs)
        print tasklist
        print len(tasklist)

    else:
        tasklist = None
        chunks = None

    tasklist = comm.bcast(tasklist, root=0)
    chunks   = comm.bcast(chunks, root=0)

    nTasks  = 1.0*len(tasklist)
    nSteps = np.ceil(nTasks/nProcs).astype(int)


    if (rank < len(chunks)):
        for ray_datenum in chunks[rank]:
            
            for plane in planes:
                # plane = row[1]

                Pdyn, ByIMF, BzIMF, W = load_TS_params(ray_datenum)
                Dst = load_Dst(ray_datenum)
                # Load Kp
                Kp = Kp_at(ray_datenum)
                # Load Ae
                AE = Ae_at(ray_datenum)

                yearday = '%d%03d'%(ray_datenum.year, ray_datenum.timetuple().tm_yday)
                milliseconds_day = (ray_datenum.second + ray_datenum.minute*60 + ray_datenum.hour*60*60)*1e3 + ray_datenum.microsecond*1e-3


                print "process %d doing %s"%(rank, plane)
                
                if plane=='XZ':
                    minx = -maxD*R_E*1e3
                    maxx = maxD*R_E*1e3
                    miny = 0
                    maxy = 0
                    minz = -maxD*R_E*1e3
                    maxz = maxD*R_E*1e3
                    nx = 200
                    ny = 1
                    nz = 200
                if plane=='XY':
                    minx = -maxD*R_E*1e3
                    maxx = maxD*R_E*1e3
                    miny = -maxD*R_E*1e3
                    maxy = maxD*R_E*1e3
                    minz = 0
                    maxz = 0
                    nx = 200
                    ny = 200
                    nz = 1
                if plane=='YZ':    
                    minx = 0
                    maxx = 0
                    miny = -maxD*R_E*1e3
                    maxy = maxD*R_E*1e3
                    minz = -maxD*R_E*1e3
                    maxz = maxD*R_E*1e3
                    nx = 1
                    ny = 200
                    nz = 200


                model_outfile='model_dump_mode%d_%d_%s_%d_%s.dat'%(modelnum, use_IGRF, yearday, milliseconds_day, plane)
                print model_outfile
                cmd = '%s '%os.path.join(ray_bin_dir, 'dumpmodel') +\
                        ' --modelnum=%d --yearday=%s --milliseconds_day=%d '%(modelnum, yearday, milliseconds_day) + \
                        '--minx=%g --maxx=%g '%(minx, maxx) +\
                        '--miny=%g --maxy=%g '%(miny, maxy) +\
                        '--minz=%g --maxz=%g '%(minz, maxz) +\
                        '--nx=%g --ny=%g --nz=%g '%(nx, ny, nz) +\
                        '--filename=%s '%(model_outfile) +\
                        '--ngo_configfile=%s '%('newray.in') +\
                        '--use_igrf=%g --use_tsyganenko=%g '%(use_IGRF,0) +\
                        '--tsyganenko_Pdyn=%g '%(Pdyn) +\
                        '--tsyganenko_Dst=%g '%(Dst) +\
                        '--tsyganenko_ByIMF=%g '%(ByIMF) +\
                        '--tsyganenko_BzIMF=%g '%(BzIMF) +\
                        '--tsyganenko_W1=%g '%(W[0]) +\
                        '--tsyganenko_W2=%g '%(W[1]) +\
                        '--tsyganenko_W3=%g '%(W[2]) +\
                        '--tsyganenko_W4=%g '%(W[3]) +\
                        '--tsyganenko_W5=%g '%(W[4]) +\
                        '--tsyganenko_W6=%g '%(W[5]) +\
                        '--gcpm_kp=%g'%(Kp)

                if modelnum==4:
                    cmd += '--interp_interpfile=%s '%(interpfile) +\
                        '--scattered_interp_window_scale=%d '%(scattered_interp_window_scale) +\
                        '--scattered_interp_order=%d '%(scattered_interp_order) +\
                        '--scattered_interp_exact=%d '%(scattered_interp_exact) +\
                        '--scattered_interp_local_window_scale=%d '%(scattered_interp_local_window_scale)

                # print cmd

                os.system(cmd)
                # print 'cp %s %s'%(model_outfile, os.path.join(ray_out_dir, model_outfile))
                os.system('cp %s %s'%(model_outfile, os.path.join(ray_out_dir, model_outfile)))

            # have all 3 models on the current node -- generate plot
            xlims = [-maxD, maxD]
            ylims = [-maxD, maxD]


            # Get direction to sun (GSE system - x axis points to sun)
            x_in = [1, 0, 0]
            sun = xf.gse2sm(x_in, ray_datenum)

            # --------------- Latex Plot Beautification --------------------------
            fig_width = 12  # width in inches
            fig_height = 4      # height in inches
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

            fig, ax = plt.subplots(1,3)

            p = []
            for ind, plane in enumerate(planes):

                earth = plt.Circle((0,0),1,color='0.5',alpha=1, zorder=100)
                iono  = plt.Circle((0,0),(R_E + H_IONO)/R_E, color='w',alpha=0.8, zorder=99)
                ax[ind].add_patch(earth)   
                ax[ind].add_patch(iono)

                model_outfile='model_dump_mode%d_%d_%s_%d_%s.dat'%(modelnum, use_IGRF, yearday, milliseconds_day, plane)
                d = readdump(model_outfile)

                Ne = d['Ns'][0,:,:,:].squeeze().T*1e-6
                Ne[np.isnan(Ne)] = 0
                
                # Plot direction to the sun
                if plane == 'XY':
                    sv = [sun[0], sun[1]]
                if plane == 'XZ':
                    sv = [sun[0], sun[2]]
                if plane == 'YZ':
                    sv = [sun[1], sun[2]]
                ax[ind].plot([0, sv[0]], [0, sv[1]],'w', linewidth=2, zorder=101)
                



                px = np.linspace(-maxD, maxD, 200)
                py = np.linspace(-maxD, maxD, 200)

                p.append(ax[ind].pcolormesh(px, py, np.log10(Ne), zorder=98))
                p[ind].set_clim(clims)
                
                
                # ax[ind].set_aspect('equal')
                ax[ind].set_xlim(xlims)
                ax[ind].set_ylim(ylims)
                ax[ind].set_title('%s'%plane)
                
            divider = make_axes_locatable(ax[2])
            cax = divider.append_axes("right",size="4%",pad=0.15)
            cb = plt.colorbar(p[2], cax=cax)
            cb.set_label('Electron Density (#/cm$^3$)')
            cticks = np.arange(clims[0],clims[1] + 1)
            cb.set_ticks(cticks)
            cticklabels = ['$10^{%d}$'%k for k in cticks]
            cb.set_ticklabels(cticklabels)

            ax[1].set_xlabel('L (R$_E$)')
            ax[0].set_ylabel('L (R$_E$)')
            ax[1].set_yticks([])
            ax[2].set_yticks([])

            figname='plasmasphere_mode%d_%d_%s_%d.png'%(modelnum, use_IGRF, yearday, milliseconds_day)

            fig.suptitle('Mode: %d %s Kp: %d - Dst: %d\n '%(modelnum, ray_datenum, Kp, Dst))

            fig.tight_layout()
            fig.savefig(figname, ldpi=300)
            os.system('mv %s %s'%(figname, os.path.join(ray_out_dir, figname)))

            plt.close('all')

            # clean up
            for plane in planes:
                model_outfile='model_dump_mode%d_%d_%s_%d_%s.dat'%(modelnum, use_IGRF, yearday, milliseconds_day, plane)
                os.system('rm %s'%model_outfile)





comm.Barrier()

if rank==0:
    print "-------- Finished with raytracing ---------"


