# Start the grid builder (random)

import numpy as np
import os
import commands
import subprocess
import shutil


raytracer_root = '/shared/users/asousa/software/raytracer_v1.17'
output_dir = '/shared/users/asousa/WIPP/lightning_power_study/gcpm_models/'
# print "using raytracer at %s"%raytracer_root
R_E = 6371e3

# Input params
L_max = 10
mlt = 1.0
kp = 4.0

minx = 9*R_E
maxx =  L_max*R_E
miny = 9*R_E
maxy =  L_max*R_E
minz = 9*R_E
maxz =  L_max*R_E

n_zero_altitude = 0
n_iri_pad = 0
n_initial_radial = 1000
n_initial_uniform = 200000
adaptive_nmax = 600000

yearday = 2001001
milliseconds_day = 0

initial_tol = 0.1

filename='gcpm_kp%d_%d_L%d_random_%d_%d_%d_%d_%d.txt'%(kp, yearday, L_max,n_zero_altitude, n_iri_pad, n_initial_radial, n_initial_uniform, adaptive_nmax)
filename = os.path.join(output_dir, filename)

cmd = '%s/bin/gcpm_dens_model_buildgrid_random --filename=%s '%(raytracer_root, filename) + \
      '--yearday=%d --milliseconds_day=%d --gcpm_kp=%d '%(yearday, milliseconds_day, kp) + \
      '--minx=%d --maxx=%d --miny=%d --maxy=%d --minz=%d --maxz=%d '%(minx, maxx, miny, maxy, minz, maxz) + \
      '--n_zero_altitude=%g --n_iri_pad=%g --n_initial_radial=%g '%(n_zero_altitude, n_iri_pad, n_initial_radial) +\
      '--n_initial_uniform=%g --adaptive_nmax=%g '%(n_initial_uniform, adaptive_nmax) + \
      '--initial_tol=%g '%initial_tol


print "cmd: ", cmd
os.system(cmd)