#include <lightning_power.h>
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"


int main(int argc, char *argv[]) 
{

    //                    alt, lat, lon (geomagnetic)
    double flash_pos[3] = {1, 50, 84};
    double flash_I0 = -10000;

    string ray_inp_dir = "/shared/users/asousa/WIPP/lightning_power_study/rays/globe_singleflash";
    string outfile_name = "/shared/users/asousa/WIPP/lightning_power_study/test_dump.dat";


    double f1 = 200;
    double f2 = 230;
    int iyr = 2001;
    int idoy = 001;
    int isec = 0;

    int num_freqs = floor((f2 - f1)/FREQ_STEP_SIZE);

    // Storage space for the first sweep
    vector<Vector3d> interp_points[num_freqs];
    vector<double>   interp_data[num_freqs];

    double in_lat, in_lon, avg_distance_from_flash;
    map <int, rayF> raylist;
    rayF cur_rays[8];
    double latmin, latmax, lonmin, lonmax;
    double wmin, wmax, tmax;
    double dlat, dlon, dw;
    double lat0, lon0, rad0;
    double tmp_coords[3], flash_pos_sm[3];
    time_t run_tstart, run_tend;
    int num_freqs_fine;
    double frame_area, initial_area;

    Vector3d cell_pos;

    Vector3d corners[8];
    double damping_at_corners[8];
    Vector3d centerpoint;
    double damping_avg;
    double bounding_sphere_radius;

    // Array of pointers to single-timestep frames
    rayT cur_frames[8];
    rayT prev_frames[8];


    int itime_in[2];
    int yearday = iyr*1000 + idoy;
    itime_in[0] = yearday;
    itime_in[1] = isec*1e3;

    const int NX = round((XMAX - XMIN)/GRID_STEP_SIZE);
    const int NY = round((YMAX - YMIN)/GRID_STEP_SIZE);
    const int NZ = round((ZMAX - ZMIN)/GRID_STEP_SIZE);


    double xaxis[NX];
    double yaxis[NY];
    double zaxis[NZ];

    for (int i=0; i<NX; ++i) { xaxis[i] = XMIN + i*GRID_STEP_SIZE; }
    for (int i=0; i<NY; ++i) { yaxis[i] = YMIN + i*GRID_STEP_SIZE; }
    for (int i=0; i<NZ; ++i) { zaxis[i] = ZMIN + i*GRID_STEP_SIZE; }

    // Set up output grid
    cout << "NX: " << NX << " NY: " << NY << " NZ: " << NZ << endl;
    double out_grid[NX][NY][NZ];

    for (int i=0; i<NX; ++i) { for (int j=0; j<NY; ++j) { for (int k=0; k<NZ; ++k) {
        out_grid[i][j][k] = 0;
    }}}



// --------------- Get flash input coordinates -----------------------    
    lat0 = D2R*flash_pos[1];
    lon0 = D2R*flash_pos[2];
    rad0 = flash_pos[0];
    // tmp_coords = {rad0, lat0, lon0};
    // sphcar(flash_pos, tmp_coords);
    // Get flash position in SM coords:
    pol_to_cart_d_(&lat0, &lon0, &rad0, tmp_coords);
    mag_to_sm_d_(itime_in, tmp_coords, flash_pos_sm);



 
    ostringstream tmp;
    tmp << ray_inp_dir << "/f_" << f1;
    vector <vector<double> > available_rays;
    vector <vector<double> > adjacent_rays;


    get_available_rays(tmp.str(), &available_rays);
    cout << "found " << available_rays.size() << " available rays\n";

    adjacent_rays = find_adjacent_rays(available_rays);
    cout << "found " << adjacent_rays.size() << " sets of adjacent rays\n";






    // -------------- Iterate through adjacent ray sets: -----------------
    for (vector < vector<double> >::iterator adj_row=adjacent_rays.begin(); adj_row!=adjacent_rays.end(); ++adj_row) {
        avg_distance_from_flash = 0;

        for (int i=0; i<4; ++i) {
            in_lat = (*adj_row)[2*i  ];
            in_lon = (*adj_row)[2*i+1];
            avg_distance_from_flash += haversine_distance(in_lat, in_lon, flash_pos[1], flash_pos[2]);
        }
        avg_distance_from_flash /= 4000.0;  // Average dist of the 4 corner rays, in km
                                            // (Assuming both frequencies have same ray spacing) 

        // Check if these rays are close enough to care about them:
        if ( (avg_distance_from_flash <= MAX_GROUND_DISTANCE) ) {
            // cout << "avg dist: " << avg_distance_from_flash << endl;
            // print_vector(*adj_row);

            // -------------- Load current rays ----------------------------------:
            for (int jj=0; jj<8; ++jj) {
                ostringstream cur_file;

                double f   = (jj < 4 ? f1 : f2);
                double lat = (jj < 4 ? (*adj_row)[2*jj]   : (*adj_row)[2*(jj - 4)]);
                double lon = (jj < 4 ? (*adj_row)[2*jj+1] : (*adj_row)[2*(jj - 4) + 1]);
                cur_file << ray_inp_dir   << "/f_" << f   << "/lon_"    << lon
                         << "/ray_" << f  << "_"   << lat << "_" << lon << ".ray";

                raylist = read_rayfile(cur_file.str());
                cur_rays[jj] = raylist.at(1); // First ray should have ray number "1"
         
                cur_rays[jj].in_lat = lat; cur_rays[jj].in_lon = lon;
                // calc_stix_parameters(&(cur_rays[jj]));

                // Also load the damping file: 
                cur_file.str(""); cur_file.clear();
                cur_file << ray_inp_dir   << "/f_" << f   << "/lon_"    << lon
                         << "/damp_" << f  << "_"   << lat << "_" << lon << ".ray";
                raylist = read_dampfile(cur_file.str());
                cur_rays[jj].damping = raylist.at(1).damping;

            }
            
            // if (DEBUG) {check_memory_usage();}
        
            // Find minimum and maximum frequencies, start lats, and start lons:
            wmin   = cur_rays[0].w;         wmax = cur_rays[0].w;
            latmin = cur_rays[0].in_lat;  latmax = cur_rays[0].in_lat;
            lonmin = cur_rays[0].in_lon;  lonmax = cur_rays[0].in_lon;
            tmax   = cur_rays[0].time.back();
            double in_lat, in_lon;

            for (int i=1; i < 8; i++) {
                in_lon = cur_rays[i].in_lon;
                in_lat = cur_rays[i].in_lat;
                
                if (in_lon >= 360)                  { in_lon -= 360; }
                if (cur_rays[i].w < wmin)          { wmin = cur_rays[i].w; } 
                if (cur_rays[i].w > wmax )         { wmax = cur_rays[i].w; }
                if (cur_rays[i].in_lat < latmin )  { latmin = cur_rays[i].in_lat; }
                if (cur_rays[i].in_lat > latmax )  { latmax = cur_rays[i].in_lat; }
                if (in_lon < lonmin )               { lonmin = in_lon; }
                if (in_lon > lonmax )               { lonmax = in_lon; }
                if (cur_rays[i].time.back() < tmax){ tmax = cur_rays[i].time.back(); }
            }

            tmax = min(tmax, TIME_MAX);

            // starting separation in lat, lon directions (meters)
            dlat = D2R*(R_E + H_IONO)*(latmax - latmin);
            dlon = D2R*(R_E + H_IONO)*(lonmax - lonmin)*cos(D2R*(latmax + latmin)/2.);
            dw   = wmax - wmin;

            cout << "\n----- current rays: -----\n";
            cout << "lon: " << lonmax << ", " << lonmin << " deg\n";
            cout << "lat: " << latmax << ", " << latmin << " deg\n";
            cout << "f: " << wmax/(2*PI) << ", " << wmin/(2*PI) << " Hz\n"; 
            cout << "Avg distance from flash: " << avg_distance_from_flash << " km\n";

            num_freqs_fine = max(1, (int)floor( (wmax - wmin)/(2.*PI*FREQ_STEP_SIZE )));


            // Scale the input power by dlat, dlon, dw:
            // (ray spacing may not be consistent)
            double inp_pwr = 0;

            inp_pwr = total_input_power(flash_pos_sm, flash_I0, 
                                        latmin, latmax, lonmin, lonmax, wmin, wmax, itime_in);
            cout << "input power: " << inp_pwr << " Watts\n";



            // --------------------- Interpolate + sum over output grid ----------------
            //                                ( The main event)            
            // -------------------------------------------------------------------------
            time(&run_tstart);
            // Interpolate the first frames:
            for (int zz=0; zz<8; zz++) { 
                interp_rayF(&cur_rays[zz], &(prev_frames[zz]), 0); 
                corners[zz] = prev_frames[zz].pos;
            }




            initial_area = polygon_frame_area(corners);

            // Step forward in time:
            for (double tt=TIME_STEP; tt < tmax; tt+=TIME_STEP) {
                // interpolate current frames:
                for (int zz=0; zz<8; zz++) { interp_rayF(&cur_rays[zz], &(cur_frames[zz]), tt);}




                // Find min and max values to search thru:
                // double cmin[3] = {1000};
                // double cmax[3] = {-1000};
                // for (int zz=0; zz<8; ++zz) {
                //     for (int bb=0; bb<3; ++bb) {
                //         if (cur_frames[zz].pos[bb] < cmin[bb]) {cmin[bb] = cur_frames[zz].pos[bb];}
                //         if (cur_frames[zz].pos[bb] > cmax[bb]) {cmax[bb] = cur_frames[zz].pos[bb];}
                //         if (prev_frames[zz].pos[bb] < cmin[bb]) {cmin[bb] = prev_frames[zz].pos[bb];}
                //         if (prev_frames[zz].pos[bb] > cmax[bb]) {cmax[bb] = prev_frames[zz].pos[bb];}
                //     }
                // }


                // int xmin_ind = nearest(xaxis, NX, cmin[0], false);
                // int xmax_ind = nearest(xaxis, NX, cmax[0], false);
                // int ymin_ind = nearest(yaxis, NX, cmin[1], false);
                // int ymax_ind = nearest(yaxis, NX, cmax[1], false);
                // int zmin_ind = nearest(zaxis, NX, cmin[2], false);
                // int zmax_ind = nearest(zaxis, NX, cmax[2], false);

                // int xmin_ind = floor((cmin[0] - XMIN)/GRID_STEP_SIZE);
                // int xmax_ind = int(ceil((cmax[0] - cmin[0])/GRID_STEP_SIZE + xmin_ind));
                // int ymin_ind = floor((cmin[1] - YMIN)/GRID_STEP_SIZE);
                // int ymax_ind = int(ceil((cmax[1] - cmin[1])/GRID_STEP_SIZE + ymin_ind));
                // int zmin_ind = floor((cmin[2] - ZMIN)/GRID_STEP_SIZE);
                // int zmax_ind = int(ceil((cmax[2] - cmin[2])/GRID_STEP_SIZE + zmin_ind));
                // print_array(cmin, 3);
                // print_array(cmax, 3);
                
                // cout << xmin_ind << " " << xmax_ind << " ";
                // cout << ymin_ind << " " << ymax_ind << " ";
                // cout << zmin_ind << " " << zmax_ind << " " << endl;


                int x_ind, y_ind, z_ind;
                // (Nested-loop hell)
                // Loop over frequencies:
                double ii, jj;
                double kk;
                // for (double kk=0; kk < 1; kk += 1./num_freqs_fine) {
                for (int k_ind = 0; k_ind < num_freqs; ++k_ind) {
                    kk = k_ind/num_freqs;

                    // interpolate corners over frequency axis:
                    for (int i=0; i<4; ++i) {
                        corners[i]   = cur_frames[i].pos*kk + cur_frames[i+4].pos*(1-kk);
                        corners[i+4] = prev_frames[i].pos*kk + prev_frames[i+4].pos*(1-kk);
                        damping_at_corners[i] = cur_frames[i].damping*kk + cur_frames[i+4].damping*(1-kk);
                        damping_at_corners[i+4] = prev_frames[i].damping*kk + prev_frames[i+4].damping*(1-kk);
                        // cout << "corner: " << corners[i].transpose() << endl;
                    }

                    // Get frame area for geometric factor:
                    frame_area = polygon_frame_area(corners);


                    // Just grabbing the value at the center point: Average from all corners.
                    centerpoint[0] = 0; centerpoint[1] = 0; centerpoint[2] = 0;
                    damping_avg = 0;
                    for (int i=0; i<8; ++i) {
                        centerpoint += corners[i];
                        damping_avg += damping_at_corners[i];
                    }

                    centerpoint /=8.0;
                    damping_avg /=8.0;


                    // cout << tt << endl;
                    interp_points[k_ind].push_back(centerpoint);
                    double cur_pwr = (inp_pwr/frame_area)*damping_avg*FREQ_STEP_SIZE;
                    interp_data[k_ind].push_back(cur_pwr);

                    // Quantize and add to grid
                    x_ind = nearest(xaxis, NX, centerpoint[0], false);
                    y_ind = nearest(yaxis, NX, centerpoint[1], false);
                    z_ind = nearest(zaxis, NX, centerpoint[2], false);

                    out_grid[x_ind][y_ind][z_ind] += cur_pwr;
                }
            // Step forward one frame:
            for (int zz=0; zz<8; zz++) { prev_frames[zz] = cur_frames[zz]; }        }

            }
        }

        // cout << "Starting alglib stuff..." << endl;
        // cout << interp_points[0].size() << " entries"<<endl;


        // // Next -- interpolate across the whole danged output grid:
        // alglib::rbfmodel model;
        // alglib::rbfcreate(3, 1, model);

        // // double rbf_inp_coords[interp_points[0].size()][4];
        // alglib::real_2d_array rbf_inps;
        // rbf_inps.setlength(interp_points[0].size(), 4);
        // double cmin[3] = {1000};
        // double cmax[3] = {-1000};
        // for (int i=0; i<interp_points[0].size(); ++i) {
        //     for (int j=0; j<3; ++j) {
        //         rbf_inps[i][j] = interp_points[0][i].data()[j];

        //         if (rbf_inps[i][j] > cmax[j]) {cmax[j] = rbf_inps[i][j]; }
        //         if (rbf_inps[i][j] < cmin[j]) {cmin[j] = rbf_inps[i][j]; }

        //     }
        //     rbf_inps[i][3] = interp_data[0][i];
        //     // cout << interp_points[0][i].transpose() << endl;
        // }


        // int xmin_ind = nearest(xaxis, NX, cmin[0], false);
        // int xmax_ind = nearest(xaxis, NX, cmax[0], false);
        // int ymin_ind = nearest(yaxis, NX, cmin[1], false);
        // int ymax_ind = nearest(yaxis, NX, cmax[1], false);
        // int zmin_ind = nearest(zaxis, NX, cmin[2], false);
        // int zmax_ind = nearest(zaxis, NX, cmax[2], false);


        // alglib::rbfreport rep;
        // alglib::rbfsetpoints(model, rbf_inps);
        // alglib::rbfsetalgoqnn(model, 0.5, 1);
        // // rbfsetalgomultilayer(model, 0.05, 10);
        // alglib::rbfbuildmodel(model, rep);

        // double xval, yval, zval;
        // for (int x_ind = xmin_ind; x_ind < xmax_ind; ++x_ind){
        //     for (int y_ind = ymin_ind; y_ind < ymax_ind; ++y_ind){
        //         for (int z_ind = zmin_ind; z_ind < zmax_ind; ++z_ind){
        //             xval = XMIN + (x_ind + 0.5)*GRID_STEP_SIZE;
        //             yval = YMIN + (y_ind + 0.5)*GRID_STEP_SIZE;
        //             zval = ZMIN + (z_ind + 0.5)*GRID_STEP_SIZE;

        //             cout << "interp at: " << xval << " " << yval << " " << zval << endl;

        //             out_grid[x_ind][y_ind][z_ind] = alglib::rbfcalc3(model, xval, yval, zval);

        //         }
        //     }
        // }        
        

    // Write output file:
    FILE* outfile;
    outfile = fopen(outfile_name.c_str(),"wb");
    if (outfile==NULL) {
        cout << "failed to open output file ;~;" << endl;
    } else {
        fwrite(out_grid, NX*NY*NZ*sizeof(double), 1, outfile);
    }
    fclose(outfile);

}