import matplotlib
matplotlib.use('agg')

from mpi4py import MPI
import numpy as np
import pickle
import gzip
import matplotlib.pyplot as plt

import commands
import subprocess
from random import shuffle
import os
import itertools
import random
import os
import time

import logging
import math

from ray_interp_Lshell import interp_ray_power, plot_LT, plot_lon
from partition import partition

# ------------ Start MPI -------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")
if "cluster" in host:
    host_num = int(host[-3:])
elif "nansen" in host:
    host_num = 0

fmt_str = '[%(levelname)s: ' + '%d/'%rank + host + '] %(message)s'
logging.basicConfig(level=logging.INFO,
                    format=fmt_str)


nProcs = 1.0*comm.Get_size()

# -------------- Simulation params ---------------------
ray_base ='/shared/users/asousa/WIPP/rays/2d/dayside/mode6/' 
out_base ='/shared/users/asousa/WIPP/lightning_power_study/outputs/dayside/mode6_v3/'
# subdirs   = ['kp0','kp2','kp4','kp6','kp8']
subdirs   = ['kp0','kp2','kp4','kp6','kp8']
# ray_dir = '/shared/users/asousa/WIPP/rays/2d/nightside/gcpm_kp0_flat'
tmax = 20
dt   = 0.1
flash_lats = np.arange(15,56,1)
# flash_lats = [47]
Llims = [1.2, 8]
L_step = 0.05

d_lon = 0.25
num_lons = 81

f_low = 200 #1120 #960 #600 # 200
f_hi  = 30000 #1310# 1120 #960 #1310 #600 #2500 #30000

max_dist = 1200
n_sub_freqs = 50

# mlt = 0 # nightside
mlt = 12 # dayside

# out_dir = '/shared/users/asousa/WIPP/lightning_power_study/outputs/nightside/gcpm_kp0'
# fig_dir = os.path.join(out_dir, 'figures')
# dat_dir = os.path.join(out_dir, 'data')

# -------------- Set up output directory --------------

if rank==0:
    for s in subdirs:
        out_dir = os.path.join(out_base, s)
        fig_dir = os.path.join(out_dir, 'figures')
        dat_dir = os.path.join(out_dir, 'data')
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        if not os.path.exists(fig_dir):
            os.mkdir(fig_dir)
        if not os.path.exists(dat_dir):
            os.mkdir(dat_dir)


# -------- Partition tasks for MPI --------------------
if rank == 0:
    lats = np.array(flash_lats, dtype=float)
    lats = lats.flatten()

    tasklist = []
    for s in subdirs:
        for l in lats:
            tasklist.append( (l, s))
    # tasklist = zip(lats)

    print "%d total tasks"%len(tasklist)
    shuffle(tasklist)
    chunks = partition(tasklist, nProcs)

else:
    tasklist = None
    chunks = None

tasklist = comm.bcast(tasklist, root=0)
chunks   = comm.bcast(chunks, root=0)

nTasks  = 1.0*len(tasklist)
nSteps = np.ceil(nTasks/nProcs).astype(int)


comm.Barrier()

# -------- Do each of the many things ---------------
if (rank < len(chunks)):
    for inp in chunks[rank]:
        lat = inp[0]
        subdir = inp[1]
        # print lat, subdir
        ray_dir = os.path.join(ray_base, subdir)
        out_dir = os.path.join(out_base, subdir)
        fig_dir = os.path.join(out_dir, 'figures')
        dat_dir = os.path.join(out_dir, 'data')

        # print ray_dir
        # print out_dir
        data_filename = os.path.join(dat_dir,'data_%d.pklz'%lat)
        cookie_cutter_filename = os.path.join(dat_dir,'cookie_%d.pklz'%lat)
        time_filename = os.path.join(fig_dir,'timeseries_%d.png'%lat)
        lon_filename  = os.path.join(fig_dir,'longitude_%d.png'%lat)
	
	if os.path.exists(data_filename):
		print data_filename + " exists!"
	else:
		print "doing ", data_filename
		data = interp_ray_power(ray_dir=ray_dir,
					tmax = tmax,
					flash_lat=lat,
					mlt=mlt,
					dt=dt,
					f_low=f_low,
					f_hi=f_hi,
					max_dist=max_dist,
					n_sub_freqs=n_sub_freqs,
					Llims=Llims,
					L_step=L_step,
					d_lon = d_lon,
					num_lons = num_lons
					)

	    #   data['data'] has dimensions ([n_fieldlines, n_freq_pairs, n_longitudes, n_times])

		# If you're doing frequencies, these files are huge (several gigs)
		# with open(data_filename,'w') as f:
		#     pickle.dump(data,f)

		# with gzip.open(data_filename,'wb') as f:
		#     pickle.dump(data,f)
		


		# # Make the cookie cutter:
		# d = np.sum(data['data'], axis=3)
		# d_new = dict()
		# d_new['data'] = d
		# d_new['Lshells'] = data['Lshells']
		# d_new['lons'] = data['lons']
		# d_new['freq_pairs'] = data['freq_pairs']

		# with gzip.open(cookie_cutter_filename,'wb') as f:
		#     pickle.dump(d_new,f)

		# (fieldline x frequency x longitude)
		data['spectrum'] = np.sum(data['data'],axis=3)

		# Sum over all frequencies (for plots)
		# (fieldline x longitude x time)
		data['data'] = np.sum(data['data'],axis=1)

		# Save the reduced data
		with gzip.open(data_filename,'wb') as f:
		    pickle.dump(data,f, protocol=pickle.HIGHEST_PROTOCOL)


		# Plot L-shell vs time (all frequencies)
		fig, ax = plot_LT(data)
		fig.suptitle('%d$^o$ latitude'%lat)
		fig.savefig(time_filename,ldpi=300)
		plt.close(fig)

		# Plot L-shell vs Longitude (all frequencies)
		fig, ax = plot_lon(data)
		fig.suptitle('%d$^o$ latitude'%lat)
		fig.savefig(lon_filename,ldpi=300)
		plt.close(fig)
