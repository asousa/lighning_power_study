// wipp.h
#ifndef ltp_H
#define ltp_H
#include <Eigen/Core>
#include <Eigen/Dense>  // Cross product lives here



#include <algorithm>    // std::next_permutation, std::sort

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>
#include <map>
// #include <ctime>
#include <getopt.h>
#include <consts.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <ftw.h>

// alglib includes
// #include <interpolation.h>
// #include "stdafx.h"

// Delauncy linterp includes:
// #include "delaunay_d_interp.h"

using namespace std;
using namespace Eigen;

// Structure for holding an entire rayfile (yuge)
typedef struct rayF {
    int size;                         // Number of timesteps
    double w;                         // frequency (angular)
    double stopcond;                  // Stop condition
    double nspec;                     // Number of species in plasmasphere model
    vector <double> time;             // Group time

    vector <vector <double> > pos;
    vector <vector <double> > vprel;
    vector <vector <double> > vgrel;
    vector <vector <double> > n;
    vector <vector <double> > B0;

    // Simulation time
    int iyr;
    int idoy;
    int isec;

    // Origin coordinates
    double in_radius, in_lat, in_lon, in_w;

    // Variable-length stuff (depending on number of constituents in model)
    vector <double> qs;              // species charge
    vector <double> ms;              // species mass
    vector <vector <double> > Ns;    // number density of species (m^-3)
    vector <vector <double> > nus;   // collision frequencies
    vector <double> damping;         // Damping vector (normalized to 1)
    double inp_pwr;                  // input power of ray

    // Stix parameters
    vector <double> stixP;
    vector <double> stixR;
    vector <double> stixL;
    // vector <double> stixS;
    // vector <double> stixD;
    // vector <double> stixA;
    // vector <double> stixB;    

} rayF;

// Structure for holding a single timestep of a ray (smol)
typedef struct rayT {
    double w;               // frequency (angular)
    double nspec;           // Number of species in plasmasphere model

    double time;            // Group time

    double dt;              // time and frequency bin size
    double dw;
    double dlat;
    double dlon;
    double ds;              // Length of ray segment (meters)
    Eigen::Vector3d pos;
    Eigen::Vector3d n;
    Eigen::Vector3d B0;
    
    // Variable-length stuff (depending on number of constituents in model)
    vector <double> qs;    // species charge
    vector <double> ms;    // species mass
    vector <double> Ns;    // number density of species (m^-3)
    vector <double> nus;   // collision frequencies
    
    double damping;        // Damping vector (normalized to 1)
    double inp_pwr;        // input power of ray

    // Stix parameters
    double stixP;
    double stixR;
    double stixL;
    // double stixS;
    // double stixD;
    double in_lat;
    double in_lon;

    int num_rays;           // Number of rays summed together here (for averaging)

} rayT;

// Store each cell in the spectrogram in cellT
typedef struct cellT {
  double        Lsh;        // L shell (Earth radii)
  double        lat;        // Latitude (degrees)
  double        t;          // Time (sec)
  double        f;          // Frequency (hz)
  double        pwr;        // Total power within cell
  double        damping;    // Attenuation (relative to 1)
  double        psi;        
  double        mu;
  double        stixP;
  double        stixR;
  double        stixL;
  Eigen::Vector3d      pos;
  double         num_rays;   
} cellT;




// Structure for holding parameters for a single EA segment 
// (planes perpendicular to field line of interest, where we
// will calculate scattering at)
typedef struct EA_segment {

    Eigen::Vector3d ea_norm;           // Vector normal to crossing plane
    Eigen::Vector3d ea_pos;            // Location in plane where field line intersects

    double Lsh;                        // L shell
    double radius;                     // Radius around field line to consider a crossing

    double lat;                        // Latitude (geomagnetic)
    double lon;                        // Longitude (geomagnetic)
    double dist_to_n;                  // Distance along field line to northern iono (R_E)
    double dist_to_s;                  // Distance along field line to southern iono (R_E) 

    double ftc_n;                      // flight-time constants (Walt 4.25)
    double ftc_s;
    // Precalculated stuff for scattering:
    double wh;                         // Electron gyrofrequency     (angular)
    double dwh_ds;                     // Derivative of gyrofrequency
    double alpha_lc;                   // Local loss-cone angle      (radians)
    double alpha_eq;                   // Equatorial loss-cone angle (radians)
    double ds;                         // Distance along field line between entries (m) (Walt 3.19)
    double dv_para_ds;                 // hm... good question.

    double Bo_ratio;                   // ratio of field at equator vs local
    
    double area;                       // Area of the EA segment
} EA_segment;


// Math functions
double l2_norm(vector<double> u);
double norm(double u[], int size);
vector<double> scalar_multiply(vector<double> u, double v);
double dot_product(vector<double>u, vector<double>v);
vector<double> add(vector<double>u, vector<double> v);
double signof(double val);



// Helpers
void print_vector(vector<double> u);
void print_array(double* arr, double len);

map<int, rayF> read_rayfile(string fileName);
map <int, rayF> read_dampfile(string fileName);
void write_rayfile(string fileName, map <int, rayF> raylist);


// Science!
double input_power_scaling(double* flash_loc, double* ray_loc, double mag_lat, double w, double i0);
double ionoAbsorp(float lat, long f);
float interpPt(float *xI, float *yI, int n, float xO);

void interp_ray_positions(rayT framelist[8],  double n_x, double n_y, double n_z, rayT* rayout);
void interp_ray_data(rayT framelist[8], double n_x, double n_y, double n_z, rayT* rayout);

void calc_stix_parameters(rayF* ray);
vector<EA_segment> init_EA_array(double lat, double lon, int itime_in[2], int model_number);

void dump_EA_array(vector<EA_segment> EA_array, string filename);


// ---- Coordinate transforms ----
extern "C" void geo2mag1_(int* iyr, double* xGEO, double* xMAG);

// ----- libxformd (Coordinate transform library used in the raytracer) -----
extern "C" void sm_to_geo_d_(int* itime, double* x_in, double* x_out);
extern "C" void geo_to_sm_d_(int* itime, double* x_in, double* x_out);

extern "C" void sm_to_mag_d_(int* itime, double* x_in, double* x_out);
extern "C" void mag_to_sm_d_(int* itime, double* x_in, double* x_out);

extern "C" void geo_to_mag_d_(int* itime, double* x_in, double* x_out);
extern "C" void mag_to_geo_d_(int* itime, double* x_in, double* x_out);

// (in radians)
extern "C" void cart_to_pol_d_(double* x_in, double* lat, double* lon, double* radius);
extern "C" void pol_to_cart_d_(double* lat, double* lon, double* radius, double* x_out);

extern "C" void sm_to_gsm_d_(int* itime, double* x_in, double* x_out);
extern "C" void gsm_to_sm_d_(int* itime, double* x_in, double* x_out);

// ---- My own transforms ----
// In-place cartesian / polar transforms. 
void carsph(double x[3]); 
void sphcar(double x[3]); 
void cardeg(double x[3]);
void degcar(double x[3]);

void carsph(double x[3], double x_out[3]); 
void sphcar(double x[3], double x_out[3]); 
void cardeg(double x[3], double x_out[3]); 
void degcar(double x[3], double x_out[3]); 

// In-place mapping of a data field between cartesian / polar frames.
void transform_data_sph2car(double lat, double lon, double d_in[3], double d_out[3]);
void transform_data_car2sph(double lat, double lon, double d_in[3], double d_out[3]);
void transform_data_geo2mag(int itime_in[2], double d_in[3], double d_out[3]);
void transform_data_mag2geo(int itime_in[2], double d_in[3], double d_out[3]);

double MLT(int itime[2], double lon);

bool descending_order(double a, double b);
int nearest(double arr[], int arr_len, double target, bool reverse_order);
// bool coarse_mask(rayT cur_rays[8], rayT prev_rays[8], EA_segment EA);

// bool crosses_EA(Vector3d l0, Vector3d l1, EA_segment EA_seg);

// double longitude_interval(double ra, double r0, double width_deg);

void interp_rayF(rayF* rayfile, rayT* frame, double t_target);


// vector <vector <int> > find_adjacent_rays(map <int, vector<double> > start_locs);
double haversine_distance(double latitude1, double longitude1, double latitude2, double longitude2);

vector<cellT> load_crossings(int itime_in[2], string filename);

cellT new_cell(rayT ray);
void add_cell(cellT* cell1, cellT* cell2);

// double polygon_frame_area(rayT frame[8], double weight);
double polygon_frame_area(Vector3d corners[4]);

void get_available_rays(string raypath, vector <vector<double> > *data);

vector < vector<double> > find_adjacent_rays(vector < vector<double> > start_locs);


int check_memory_usage();


double graf_iono_absorp(float lat, long f, double mlt);
double input_power_scaling(double* flash_loc, double* ray_loc,
                         double mag_lat, double w, double i0, double MLT);
double total_input_power(double flash_pos_sm[3], double i0,
                         double latmin, double latmax, 
                         double lonmin, double lonmax, 
                         double wmin, double wmax, int itime_in[2]);
 

// void find_crossing(rayT cur_frames[8], rayT prev_frames[8], Vector3d targ_point, double n_f, double *n_x, double *n_y);
void find_crossing(Vector3d corners[8], double data[8], Vector3d targ_point, double *n_x, double *n_y);

double interp_damping(rayT framelist[8], double n_x, double n_y, double n_z);

// double interp_damping_inv_dist(rayT cur_frames[8], rayT prev_frames[8], double freq_weight, Vector3d targ_point);
double interp_damping_inv_dist(Vector3d corners[8], double data_at_corners[8], Vector3d targ_point);

double bounding_sphere(Vector3d corners[8], double center_arr[3]);


#endif
