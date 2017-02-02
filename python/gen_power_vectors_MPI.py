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

# Here's how to run the c code:
# (this code interpolates rays + damping onto a uniform time axis,
# then calculates geometric spreading at the center of each "tube", defined
# by four corner rays. The result is a measure of wave power density (w/m^2) at
# each timestep. )
flower = 200; flupper = 30000;
num_freqs = 33
flogs = np.linspace(np.log10(flower), np.log10(flupper), num_freqs)
freqs = np.round(pow(10, flogs)/10.)*10
effs = zip(freqs[0:-1], freqs[1:])


ray_inp_dir = '/shared/users/asousa/WIPP/lightning_power_study/rays/nightside/ngo_igrf'
file_out_dir = '/shared/users/asousa/WIPP/lightning_power_study/outputs/nightside/ngo_igrf'
time_max = 15
num_times = 150
max_ground_distance = 2000
flash_lat = 50
flash_lon = 84

ideal_step_size = 500 # hz


if not os.path.exists(file_out_dir):
    os.mkdir(file_out_dir)

for f_ind, ff in enumerate(effs):
    
    num_freqs = int(max(1, np.round((ff[1] - ff[0])/ideal_step_size)))

    outfile_name= os.path.join(file_out_dir, 'power_vectors_%d_%d.dat'%(ff[0], ff[1]))

    cmd = '../bin/calc_power --ray_dir=%s --out_file=%s'%(ray_inp_dir, outfile_name) + \
            ' --f1=%d --f2=%d --t_max=%g --num_times=%d'%(ff[0], ff[1], time_max, num_times) + \
            ' --max_dist=%d --num_freqs=%d'%(max_ground_distance, num_freqs) + \
            ' --lat=%d --lon=%d'%(flash_lat, flash_lon)

    print cmd


    ret_val = os.system(cmd)

    if ret_val != 0:
        print "C code threw some error... Hmm."
    else:
        print "finished %d of %d"%(f_ind +1 , len(effs))