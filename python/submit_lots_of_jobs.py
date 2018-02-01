import os

# inlats = [15, 25, 35, 45, 55]
inlats = [20, 30, 40, 50]
# inlats = [40]
# inlats = [15, 25, 35, 45, 55]
kps = [0, 2, 4, 6, 8]

for inlat in inlats:
    for kp in kps:
        cmd = 'qsub -N frames_%d_%d'%(inlat, kp)  + ' -v inlat=%d'%inlat +\
            ',kp=%d jobs/movieframes.pbs'%kp
        os.system(cmd)