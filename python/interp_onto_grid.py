import numpy as np
from scipy import interpolate
from scipy.spatial import Delaunay

# for plotting
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable





def interp_onto_grid(fname, xlims, ylims, zlims, step_size, method='linear'):

    


    R_E = 6371e3
    H_IONO = 1000e3
    
    data_raw = np.loadtxt(fname)

    nf = int(data_raw[0]) # num frequencies
    nv = int(data_raw[1]) # num vectors

    # Simple file (all vectors the same length)
    # dgrid = data_raw[3:]
    # dgrid = np.reshape(dgrid, [nf,nv,nt,4],order='c')

    dlengths = np.array([int(x) for x in data_raw[2:(2 + nf*nv)]])
    nt = max(dlengths)
    dgrid = np.zeros([nf, nv, max(dlengths), 4])
    startind = 2 + nf*nv
    stopind  = 0 
    for f_ind in range(nf):
        for v_ind in range(nv):
            d_ind = f_ind*nv + v_ind
            dl = dlengths[d_ind]      
    #         print f_ind, v_ind, d_ind, dl
            stopind = startind + 4*dl
            tmp = data_raw[startind:stopind]
            curdata = np.reshape(tmp, [dl,4])
            dgrid[f_ind, v_ind, 0:dl, :] = curdata
            startind=stopind


    # print dgrid[0,:,0,0:4]


    xx = np.arange(xlims[0], xlims[1], step_size)
    yy = np.arange(ylims[0], ylims[1], step_size)
    zz = np.arange(zlims[0], zlims[1], step_size)

    nx = len(xx) 
    ny = len(yy)
    nz = len(zz)
    print nf, nv


    # Output space
    data = np.zeros([nx, ny, nz])

    starting_coords = dgrid[0,:,0,0:3]

    sc_polar = np.array([c2p3(np.array(row)) for row in starting_coords])
    # Set radius to 1, round to 1-degree bins
    sc_polar[:,0] = 1
    sc_polar[:,1:]*=180./np.pi
    sc_polar = np.round(sc_polar)

    # Delaunay triangulation on lat, lon. 
    # (This is a little sloppy since we're taking cartesian distances from
    # angular values, but fuck it, it works in this case. Might break elsewhere.)
    tris = Delaunay(sc_polar[:,1:3])
    adj_inds = np.array(tris.simplices)

    for f_ind in range(nf):
        for t in range(0,nt-1):
            hits = np.zeros([nx, ny, nz])
            for adj_row in adj_inds:
                # Check max length of each vector:
                cur_t = min(dlengths[adj_row])
                
                if t < cur_t:
                    points = dgrid[f_ind,adj_row,t:t+2,0:3]
                    vals = dgrid[f_ind,adj_row,t:t+2,3]
                    vals_flat = vals.ravel()
                    points_flat = np.vstack([points[:,:,0].ravel(), points[:,:,1].ravel(), points[:,:,2].ravel()]).T
                    
                    # Narrow down our search space to just a box around the current hull:
                    minx = min(points_flat[:,0])
                    maxx = max(points_flat[:,0])
                    miny = min(points_flat[:,1])
                    maxy = max(points_flat[:,1])
                    minz = min(points_flat[:,2])
                    maxz = max(points_flat[:,2])

                    ix = np.where((xx >= minx) & (xx < maxx))[0]
                    iy = np.where((yy >= miny) & (yy < maxy))[0]
                    iz = np.where((zz >= minz) & (zz < maxz))[0]

                    px, py, pz = np.meshgrid(ix, iy, iz, indexing='ij')  # in 3d, ij gives xyz, xy gives yxz. dumb.

                    newpoints = np.vstack([xx[px.ravel()], yy[py.ravel()], zz[pz.ravel()]]).T

                    if method=='linear':
                        # If we don't specify a fill value, 'linear' mode returns NaN for anything outside
                        # the convex hull of the current point cloud (which is ideal -- we don't want any
                        # values calculated outside the (non-convex) hull.)
                        tmp_data = interpolate.griddata(points_flat, vals_flat, newpoints, method='linear', rescale=False)
                        tmp_data = tmp_data.reshape([len(ix), len(iy), len(iz)])

                        # tmp_data = np.unravel_index
                        isnans = np.isnan(tmp_data)
                        tmp_data[np.isnan(tmp_data)] = 0
            #             print np.sum(isnans), np.sum(~isnans)
                        data[px,py,pz] += tmp_data
                        hits[px,py,pz] += ~isnans
                    elif method=='mean':
                        # fill each voxel with the mean value of all corners.
                        newtri = Delaunay(points_flat)
                        hit_mask = newtri.find_simplex(newpoints) >=0

                        hit_mask = hit_mask.reshape([len(ix), len(iy), len(iz)])
                        hits[px,py,pz] += hit_mask
                        data[px[hit_mask],py[hit_mask], pz[hit_mask]] += np.mean(vals_flat)

                # print np.max(hits), np.sum(hits > 1), len(hits)
            # Average any bins which got more than one hit at this timestep:
            # (this should just be the edges and corners)
            data[hits!=0] /= hits[hits!=0]  

    return data





# Plot average power on each axis:
def plot_avg_power_3up(data, xlims, ylims, zlims, step_size):
    R_E = 6371e3
    H_IONO = 1000e3

    xx = np.arange(xlims[0], xlims[1], step_size)
    yy = np.arange(ylims[0], ylims[1], step_size)
    zz = np.arange(zlims[0], zlims[1], step_size)

    nx = len(xx) 
    ny = len(yy)
    nz = len(zz)

    xz_sum = np.log10((np.sum(data, axis=1).T)/len(yy))
    xy_sum = np.log10((np.sum(data, axis=2).T)/len(zz))
    yz_sum = np.log10((np.sum(data, axis=0).T)/len(xx))
    
    # print "xz min, max: ", np.min(xz_sum), np.max(xz_sum)

    maxlog = np.ceil(np.max([xz_sum, xy_sum, yz_sum]))

    # Show about 3 orders of magnitude
    clims = [maxlog - 3, maxlog]
    xy_sum[np.isinf(xy_sum)] = -100
    xz_sum[np.isinf(xz_sum)] = -100
    yz_sum[np.isinf(yz_sum)] = -100



    fig, ax = plt.subplots(1,3)
    # Plot the earth
    for i in [0, 1, 2]:
        earth = plt.Circle((0,0),1,color='0.5',alpha=1, zorder=100)
        iono  = plt.Circle((0,0),(R_E + H_IONO)/R_E, color='w',alpha=0.8, zorder=99)
        ax[i].add_patch(earth)   
        ax[i].add_patch(iono)
    
    p0 = ax[0].pcolorfast(xx, yy, xy_sum, vmin=clims[0], vmax=clims[1])
    p1 = ax[1].pcolorfast(xx, zz, xz_sum, vmin=clims[0], vmax=clims[1])
    p2 = ax[2].pcolorfast(yy, zz, yz_sum, vmin=clims[0], vmax=clims[1])

    divider = make_axes_locatable(ax[2])
    cax = divider.append_axes("right",size="4%",pad=0.15)
    cb = plt.colorbar(p2, cax=cax)
    cb.set_label('avg wave power density')
    cticks = np.arange(clims[0],clims[1] + 1)
    cb.set_ticks(cticks)
    cticklabels = ['$10^{%d}$'%k for k in cticks]
    cb.set_ticklabels(cticklabels)
    
#     ax[0].set_aspect('equal')
#     ax[1].set_aspect('equal')
#     ax[2].set_aspect('equal')

    ax[0].set_title('XY')
    ax[1].set_title('XZ')
    ax[2].set_title('YZ')
    
    ax[0].set_xlim([xx[0],xx[-1]])
    ax[0].set_ylim([yy[0],yy[-1]])
    ax[1].set_xlim([xx[0],xx[-1]])
    ax[1].set_ylim([zz[0],zz[-1]])
    ax[2].set_xlim([yy[0],yy[-1]])
    ax[2].set_ylim([zz[0],zz[-1]])
    
    
    

    
    fig.tight_layout()





def c2p3(din):
    # Cartesian to spherical (rad). Matches the fortran CART_TO_POL method
    dout = np.zeros_like(din)

    row=np.sqrt(din[0]**2+din[1]**2)
    dout[0]=np.sqrt(row*row+din[2]**2)
    dout[1]=np.arctan2(din[2],row)
    dout[2]=np.arctan2(din[1],din[0])
    dout[2]+= + (1.0 -np.sign(dout[2]))*np.pi
    
    
#     dout = np.zeros_like(din)
#     dout[0] = np.linalg.norm(din)
#     dout[1] = np.arctan2(din[1],din[0])
#     dout[2] = np.arctan2(np.sqrt(din[0]**2 + din[1]**2),din[2])
    return dout

def p2c3(din):
    # Spherical to cartesian (rad)
    dout = np.zeros_like(din)
    dout[0] = din[0]*np.cos(din[1])*np.cos(din[2])
    dout[1] = din[0]*np.cos(din[1])*np.sin(din[2])
    dout[2] = din[0]*np.sin(din[1])
    return dout