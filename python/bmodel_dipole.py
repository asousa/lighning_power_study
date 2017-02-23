import numpy as np
import igrf12
from sciencedates import datetime2yeardec

def bmodel_dipole(x_in):
    Bo = (3.12e-5)
    D2R = np.pi/180.
    R_E = 6371e3;

    R = x_in[0]
    theta = (90. - x_in[1])*D2R
    phi   = x_in[2]*D2R
    
    Bor3  = Bo*pow(R/R_E, -3.0)
    
    Brad = -2.0*Bor3*np.cos(theta);
    Btheta = -1.0*Bor3*np.sin(theta);
    Bphi = 0.0;    # Dipole model has no variation in longitude (here for completeness)

    B_out = np.zeros(3)
    B_out[0] = Brad;        # Up
    B_out[1] = Btheta;      # South
    B_out[2] = Bphi;        # East

    return B_out


def bmodel_igrf(x_in, itime):
    print "herpy derpy"
    isv=0  # Field = 0, secular variaton = 1
    itype = 1 # 1=geodetic earth, 2 = geocentric earth

    yeardec = datetime2yeardec(dtime)
    colat,elon = latlon2colat(glat,glon)

    x = np.empty(colat.size);  y = np.empty_like(x); z = np.empty_like(x); f=np.empty_like(x)
    for i,(clt,eln) in enumerate(nditer((colat,elon))):
        x[i],y[i],z[i],f[i] = igrf12.igrf12syn(isv, yeardec, itype, alt, clt, eln)

    return x.reshape(colat.shape), y.reshape(colat.shape), z.reshape(colat.shape),f.reshape(colat.shape), yeardec

def latlon2colat(glat,glon):
    #atleast_1d for iteration later
    colat = 90-np.atleast_1d(glat)
    elon = (360 + np.atleast_1d(glon)) % 360
    return colat,elon
