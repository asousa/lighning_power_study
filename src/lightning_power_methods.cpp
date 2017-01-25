#include <lightning_power.h>



void interp_rayF(rayF* rayfile, rayT* frame, double t_target) {
// float interpPt(float *xI, float *yI, int n, float t_target )
                //  time axis   data vector length   t_targ       
    int i, iHigh, iLow, iMid;
    double yO;

    double M;
    vector <double> xI = rayfile->time;
    int n = rayfile->time.size();


    // Check that t_target is within bounds
    if( (t_target < xI[0]) || t_target  > xI[n-1] ) {
        printf("\nPoint is out of bounds! %g, {%g, %g}\a\n",t_target , xI[0], xI[n-1]);
        return;
    }
      
    // Do a binary search for the correct index 
    iHigh = n-1;
    iLow = 0;  
    while(1) {
        iMid = (iHigh+iLow)/2;
        if( (t_target  >= xI[iLow]) && (t_target  < xI[iMid]) ) {
            iHigh = iMid;
        } else {
            iLow = iMid;
        }
        if(t_target ==xI[n-1]){printf("\nin interpPt\n"); return;}
        if(iHigh==iLow) {
            printf("\nexiting from binary search in 1st condtion\n");
            break;
        }
        if( (t_target >=xI[iMid]) && (t_target <xI[iMid+1]) ) break;
        if( (t_target >=xI[iMid-1]) && (t_target <xI[iMid]) ) {
           iMid--;
           break;
            }
        }

    M = ( t_target -xI[iMid] )/( xI[iMid+1]-xI[iMid] );
    // Now, let's interpolate the output values:
    // Vector-valued
    for (int k = 0; k < 3; k++) {
        frame->pos[k] = ( rayfile->pos[iMid+1][k]-rayfile->pos[iMid][k] )*M + rayfile->pos[iMid][k];
        frame->n[k] =   ( rayfile->n[iMid+1][k]  -rayfile->n[iMid][k]   )*M + rayfile->n[iMid][k];
        frame->B0[k] =  ( rayfile->B0[iMid+1][k] -rayfile->B0[iMid][k]  )*M + rayfile->B0[iMid][k];
    }
    // Scalar-valued
    frame->damping = ( rayfile->damping[iMid+1]-rayfile->damping[iMid] )*M + rayfile->damping[iMid];
    // frame->stixP = ( rayfile->stixP[iMid+1]-rayfile->stixP[iMid] )*M + rayfile->stixP[iMid];
    // frame->stixR = ( rayfile->stixR[iMid+1]-rayfile->stixR[iMid] )*M + rayfile->stixR[iMid];
    // frame->stixL = ( rayfile->stixL[iMid+1]-rayfile->stixL[iMid] )*M + rayfile->stixL[iMid];
    
    // Stuff that doesn't need interpolation:
    frame->w    = rayfile->w;
    frame->in_lat = rayfile->in_lat;
    frame->in_lon = rayfile->in_lon;
    frame->time = t_target;
    frame->inp_pwr = rayfile->inp_pwr;
}



double polygon_frame_area(rayT frame[8], double weight) {
    // Calculates the area enclosed by the set of guide rays.
    // There's two sloppy hacks here:
    //      1: the area cross product depends on the order of
    //      the indices - I'm just doing all permutations, and
    //      picking the largest entry.
    //      2: we assume the lower-frequency corners are the
    //      first four, and the upper-frequency are the second.
    // 
    //      weight is the interpolating quantity, ranging 0 to 1

    Vector3d cp(0,0,0);
    int inds[4] = {0,1,2,3};
    int n=4;
    
    double max_area_low = 0;
    double max_area_high= 0;
    double area;

    // Lower frequency rays are the first four:
    do {
        area = 0;
        Vector3d cp(0,0,0);

        for (int i=0; i<n; ++i) {
            Vector3d v1 = frame[inds[i]].pos;
            Vector3d v2 = frame[inds[(i+1)%n]].pos;
            cp += v1.cross(v2);
        }

        area = pow(R_E, 2)*cp.norm()/2.;

        if (area > max_area_low) { max_area_low = area;}
    } while ( next_permutation(inds,inds + 2) );


    // upper frequency rays are the second four:
    inds = {4,5,6,7};
    do {
        area = 0;
        Vector3d cp(0,0,0);

        for (int i=0; i<n; ++i) {
            Vector3d v1 = frame[inds[i]].pos;
            Vector3d v2 = frame[inds[(i+1)%n]].pos;
            cp += v1.cross(v2);
        }

        area = pow(R_E, 2)*cp.norm()/2.;

        if (area > max_area_high) { max_area_high = area;}
    } while ( next_permutation(inds,inds + 2) );

    // cout << max_area_low << " " << max_area_high << endl;
    // cout << max_area_high/max_area_low << endl;
    return weight*max_area_low + (1.0 - weight)*max_area_high;
}



double interp_damping(rayT framelist[8], double n_x, double n_y, double n_z) {
    // Interpolate quantities at steps nx, ny, nz.
    // (just doing damping right now - I think that's all we'll need)

    double W[8];
    W[0] = (1. - n_x)*(1. - n_y)*(1. - n_z);
    W[1] = n_x*(1. - n_y)*(1. - n_z);
    W[2] = n_y*(1. - n_x)*(1. - n_z);
    W[3] = n_x*n_y*1.*(1. - n_z);
    W[4] = (1. - n_x)*(1. - n_y)*n_z;
    W[5] = n_x*(1. - n_y)*n_z;
    W[6] = n_y*(1. - n_x)*n_z;
    W[7] = n_x*n_y*n_z*1.;

    double damping = 0;

    for (int jj=0; jj<8; jj++){  // Corner rays
        damping += W[jj]*(framelist[jj].damping);
    }

    return damping;

}

