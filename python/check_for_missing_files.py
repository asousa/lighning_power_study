import numpy as np

import os

import itertools
import random
import time

import commands

import shutil
from random import shuffle


inp_lats = np.arange(10, 61, 5) #[35] #np.arange(30, 61, 5) #[40, jh41, 42, 43]
# inp_lons = np.arange(0, 360, 5)

# Nightside -- at midnight, nightside is ~84 deg longitude.
inp_lons = np.arange(65, 106, 2)

f1 = 200; f2 = 30000;
num_freqs = 33
flogs = np.linspace(np.log10(f1), np.log10(f2), num_freqs)
freqs = np.round(pow(10, flogs)/10.)*10


project_root = '/shared/users/asousa/WIPP/lightning_power_study/'
raytracer_root = '/shared/users/asousa/software/raytracer_v1.17/'
damping_root = '/shared/users/asousa/software/damping/'
ray_bin_dir    = os.path.join(raytracer_root, 'bin')
ray_out_dir = os.path.join(project_root, 'rays','nightside','ngo_dipole')

missing_rays = 0
missing_damp = 0
for freq in freqs:
    for lat in inp_lats:
        for lon in inp_lons:
            ray_out_subdir = os.path.join(ray_out_dir, "f_%d"%freq, "lon_%d"%(lon))
            ray_outfile   = os.path.join(ray_out_subdir, 'ray_%d_%d_%d.ray'%(freq, lat, lon))
            damp_outfile  = os.path.join(ray_out_subdir,'damp_%d_%d_%d.ray'%(freq, lat, lon))

            if not (os.path.exists(ray_outfile)):
                print "missing rayfile at %d, %d, %d"%(freq, lat, lon)
                missing_rays += 1
            if not (os.path.exists(damp_outfile)):
                print "missing dampfile at %d, %d, %d"%(freq, lat, lon)
                missing_damp += 1

print "Missing %d rays"%missing_rays
print "Missing %d damp"%missing_damp
