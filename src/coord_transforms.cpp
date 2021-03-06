#include <lightning_power.h>
// Various transforms between coordinate frames.
// Includes wrappers to whichever Fortran libraries we're using,
// so that the interface is a little nicer in the rest of the code.


void cardeg(double x[3]) {
// Cartesian to polar (degrees)
    carsph(x);
    x[1] = R2D*x[1];
    x[2] = R2D*x[2];
}

void cardeg(double x_in[3], double x_out[3]) {
    carsph(x_in, x_out);
    x_out[1]*= R2D;
    x_out[2]*= R2D;
}

void degcar(double x[3]) {
// Polar to Cartesian (degrees)
    x[1] = D2R*x[1];
    x[2] = D2R*x[2];
    sphcar(x);
}

void degcar(double x[3], double x_out[3]) {
    x_out[0] = x[0]*cos(x[1]*D2R)*cos(x[2]*D2R);
    x_out[1] = x[0]*cos(x[1]*D2R)*sin(x[2]*D2R);
    x_out[2] = x[0]*sin(x[1]*D2R);
}
void carsph(double x[3]) {
    // in-place rotation from Cartesian to Spherical (radians)
    // output is R, Theta, Phi

    double lat, lon, r;
    
    r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    lon = atan2(x[1], x[0]);
    lat = asin(x[2]/r);

    x[0] = r;
    x[1] = lat;
    x[2] = lon;
}

void carsph(double x[3], double x_out[3]) {
    // cartesian to spherical, not in-place.
    double lat, lon, r;
    
    x_out[0] = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    x_out[2] = atan2(x[1], x[0]);
    x_out[1] = asin(x[2]/x_out[0]);
}

void sphcar(double x[3]) {
    // in-place rotation from Spherical to Cartesian (radians)
    // input is R, Theta, Phi, output is x y z
    double lat, lon, r;
    r = x[0]; lat = x[1]; lon = x[2];
    x[0] = r*cos(lat)*cos(lon);
    x[1] = r*cos(lat)*sin(lon);
    x[2] = r*sin(lat);
}

void sphcar(double x[3], double x_out[3]) {
    x_out[0] = x[0]*cos(x[1])*cos(x[2]);
    x_out[1] = x[0]*cos(x[1])*sin(x[2]);
    x_out[2] = x[0]*sin(x[1]);
}

void transform_data_sph2car(double lat, double lon, double d_in[3], double d_out[3]) {
    // Map a vector field from spherical to cartesian coordinates
    double M[3][3];
    
    d_out[0] = 0; d_out[1] = 0; d_out[2] = 0;

    double theta = D2R*(90. - lat);
    double phi = D2R*lon;

    double st = sin(theta);
    double sp = sin(phi);
    double ct = cos(theta);
    double cp = cos(phi);

    // Transformation matrix
    M[0][0] = st*cp;    M[0][1] = ct*cp;   M[0][2] = -sp;
    M[1][0] = st*sp;    M[1][1] = ct*sp;   M[1][2] = cp;
    M[2][0] = ct;       M[2][1] = -st;     M[2][2] = 0;

    // Matrix multiply
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            d_out[row] += M[row][col]*d_in[col];
        }
    }
}

void transform_data_car2sph(double lat, double lon, double d_in[3], double d_out[3]) {
    // Map a vectyor field from cartesian to spherical coordinates
 double M[3][3];
    
    d_out[0] = 0; d_out[1] = 0; d_out[2] = 0;

    double theta = D2R*(90. - lat);
    double phi = D2R*lon;

    double st = sin(theta);
    double sp = sin(phi);
    double ct = cos(theta);
    double cp = cos(phi);

    // Transformation matrix
    M[0][0] = st*cp;    M[0][1] = st*sp;   M[0][2] = ct;
    M[1][0] = ct*cp;    M[1][1] = ct*sp;   M[1][2] = -st;
    M[2][0] = -sp;      M[2][1] = cp;      M[2][2] = 0;

    // Matrix multiply
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            d_out[row] += M[row][col]*d_in[col];
        }
    }
}

void transform_data_geo2mag(int itime_in[2], double d_in[3], double d_out[3]) {
    // Map data from geographic (cartesian) to geomagnetic (cartesian)
    // Basis vectors in input frame:
    double A1[3] = {1, 0, 0};
    double A2[3] = {0, 1, 0};
    double A3[3] = {0, 0, 1};

    // Basis vectors in output frame:
    double B1[3];
    double B2[3];
    double B3[3];

    double data_mag[3];

    geo_to_mag_d_(itime_in, A1, B1);
    geo_to_mag_d_(itime_in, A2, B2);
    geo_to_mag_d_(itime_in, A3, B3);

    // Inner product: 
    d_out[0] = d_in[0]*B1[0] + d_in[1]*B2[0] + d_in[2]*B3[0];
    d_out[1] = d_in[0]*B1[1] + d_in[1]*B2[1] + d_in[2]*B3[1];
    d_out[2] = d_in[0]*B1[2] + d_in[1]*B2[2] + d_in[2]*B3[2];
}

void transform_data_mag2geo(int itime_in[2], double d_in[3], double d_out[3]) {
    // Map data from geomagnetic (cartesian) to geographic (cartesian)
    // Basis vectors in input frame:
    double A1[3] = {1, 0, 0};
    double A2[3] = {0, 1, 0};
    double A3[3] = {0, 0, 1};

    // Basis vectors in output frame:
    double B1[3];
    double B2[3];
    double B3[3];

    double data_mag[3];

    mag_to_geo_d_(itime_in, A1, B1);
    mag_to_geo_d_(itime_in, A2, B2);
    mag_to_geo_d_(itime_in, A3, B3);

    // Inner product: 
    d_out[0] = d_in[0]*B1[0] + d_in[1]*B2[0] + d_in[2]*B3[0];
    d_out[1] = d_in[0]*B1[1] + d_in[1]*B2[1] + d_in[2]*B3[1];
    d_out[2] = d_in[0]*B1[2] + d_in[1]*B2[2] + d_in[2]*B3[2];
}


double haversine_distance(double latitude1, double longitude1, double latitude2,
                          double longitude2) {
    // Great-circle distance between two pairs of (lat, lon), in degrees.
    double lat1 = D2R*latitude1;
    double lon1 = D2R*longitude1;
    double lat2 = D2R*latitude2;
    double lon2 = D2R*longitude2;

    double d_lat = fabs(lat1 - lat2);
    double d_lon = fabs(lon1 - lon2);

    double a = pow(sin(d_lat / 2), 2) + cos(lat1) * cos(lat2) * pow(sin(d_lon / 2), 2);

    //double d_sigma = 2 * atan2(sqrt(a), sqrt(1 - a));
    double d_sigma = 2 * asin(sqrt(a));

    return R_E * d_sigma;
}


double MLT(int itime[2], double lon) {
    // Input: itime, lon in geomagnetic dipole coords.
    // Ref: "Magnetic Coordinate Systems", Laundal and Richmond
    // Space Science Review 2016, DOI 10.1007/s11214-016-0275-y 

    double ut_hr = itime[1]/1000.0/60.0;  // Milliseconds to fractional hours (UT)
    double A1[3] = {1, 51.48, 0};         // Location of Greenwich (for UT reference) 
    double B1[3];                         // Location of Greenwich in geomag

    degcar(A1);
    geo_to_mag_d_(itime, A1, B1);
    cardeg(B1);

    return fmod(ut_hr + (lon - B1[2])/15.0,  24);
}
