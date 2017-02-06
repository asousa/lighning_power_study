from mpi4py import MPI
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
import datetime as dt
import commands
from partition import partition


# Here's how to run the c code:
# (this code interpolates rays + damping onto a uniform time axis,
# then calculates geometric spreading at the center of each "tube", defined
# by four corner rays. The result is a measure of wave power density (w/m^2) at
# # each timestep. )
flower = 200; flupper = 30000;
num_freqs = 33
flogs = np.linspace(np.log10(flower), np.log10(flupper), num_freqs)
freqs = np.round(pow(10, flogs)/10.)*10

effs = zip(freqs[0:-1], freqs[1:])



ray_inp_dir = '/shared/users/asousa/WIPP/rays/nightside/ngo_igrf'
file_out_dir = '/shared/users/asousa/WIPP/lightning_power_study/outputs/testing/'
run_path = '/shared/users/asousa/WIPP/lightning_power_study/bin/calc_power'

time_max = 20
num_times = 200
max_ground_distance = 2000
# flash_lats = [20, 30, 40, 50]
flash_lats = [21, 31, 41, 51]
flash_lon = 76          # antisolar point at 1/1/2010:0:0:0 ~ 25.7 lat, 76.6 lon
ideal_step_size = 500 # hz


# ------------ Start MPI -------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")

nProcs = 1.0*comm.Get_size()

if rank == 0:
    tasklist = effs
    # tasklist = [(launch_alt, w,x,y) for w,x,y in itertools.product(inp_lats, inp_lons, freqs)]
    chunks = partition(tasklist, nProcs)

else:
    tasklist = None
    chunks = None

tasklist = comm.bcast(tasklist, root=0)
chunks   = comm.bcast(chunks, root=0)

nTasks  = 1.0*len(tasklist)
nSteps = np.ceil(nTasks/nProcs).astype(int)

if (rank < len(chunks)):
    # working_path = os.path.join(os.path.expanduser("~"),"rayTmp_%d"%(rank))
    print "Subprocess %s on %s: doing %d tasks"%(rank, host, len(chunks[rank]))
    for ff in chunks[rank]:
        for flash_lat in flash_lats:

            outfile_name= os.path.join(file_out_dir, 'lat_%d'%flash_lat, 'power_vectors_%d_%d.dat'%(ff[0], ff[1]))
            num_freqs = int(max(1, np.round((ff[1] - ff[0])/ideal_step_size)))

            cmd = '%s --ray_dir=%s --out_file=%s'%(run_path, ray_inp_dir, outfile_name) + \
                    ' --f1=%d --f2=%d --t_max=%g --num_times=%d'%(ff[0], ff[1], time_max, num_times) + \
                    ' --max_dist=%d --num_freqs=%d'%(max_ground_distance, num_freqs) + \
                    ' --lat=%g --lon=%g'%(flash_lat, flash_lon)

            print cmd
            ret_val = os.system(cmd)

            if ret_val != 0:
                print "C code threw some error... Hmm."
            else:
                print "finished entry between %d and %d"%(ff[0], ff[1])








# file_out_dir = '/shared/users/asousa/WIPP/lightning_power_study/outputs/nightside/ngo_igrf/lat_%d'%flash_lat

# if not os.path.exists(file_out_dir):
#     os.mkdir(file_out_dir)

# for f_ind, ff in enumerate(effs):
    
#     num_freqs = int(max(1, np.round((ff[1] - ff[0])/ideal_step_size)))

#     outfile_name= os.path.join(file_out_dir, 'power_vectors_%d_%d.dat'%(ff[0], ff[1]))

#     cmd = '../bin/calc_power --ray_dir=%s --out_file=%s'%(ray_inp_dir, outfile_name) + \
#             ' --f1=%d --f2=%d --t_max=%g --num_times=%d'%(ff[0], ff[1], time_max, num_times) + \
#             ' --max_dist=%d --num_freqs=%d'%(max_ground_distance, num_freqs) + \
#             ' --lat=%d --lon=%d'%(flash_lat, flash_lon)

#     print cmd


#     ret_val = os.system(cmd)

#     if ret_val != 0:
#         print "C code threw some error... Hmm."
#     else:
#         print "finished %d of %d"%(f_ind +1 , len(effs))