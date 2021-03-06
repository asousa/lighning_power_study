#include <lightning_power.h>
// #include "stdafx.h"
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include "interpolation.h"


int main(int argc, char *argv[]) 
{




    // Default arguments:
    //                    alt, lat, lon (geomagnetic)
    double flash_pos[3] = {1, 50, 84};
    double flash_I0 = -10000;

    string ray_inp_dir = "/shared/users/asousa/WIPP/lightning_power_study/rays/globe_ngo";
    string outfile_name = "/shared/users/asousa/WIPP/lightning_power_study/test_dump.dat";


    double f1 = 200;
    double f2 = 230;
    int iyr = 2010;
    int idoy = 001;
    int isec = 0;
    double time_max = 10.0;             // Sec
    int num_times = 100;                  
    double max_ground_distance = 1000;  // Km
    int num_freqs = 1;

    // double freq_step_size = 100;        // Hz
    // Parse input arguments:


    int opt;
    int opt_index;
    static struct option long_options[] =
    {
        {"out_file",      required_argument,    0, 'a'},
        {"iyr",           required_argument,    0, 'b'},
        {"idoy",          required_argument,    0, 'c'},
        {"isec",          required_argument,    0, 'd'},
        {"I0",            required_argument,    0, 'e'},
        {"f1",            required_argument,    0, 'f'},
        {"f2" ,           required_argument,    0, 'g'},
        {"ray_dir",       required_argument,    0, 'h'},
        {"t_max",         required_argument,    0, 'i'},
        {"max_dist",      required_argument,    0, 'j'},
        {"num_times",     required_argument,    0, 'k'},
        {"num_freqs",     required_argument,    0, 'l'},
        {"lat",           required_argument,    0, 'm'},
        {"lon",           required_argument,    0, 'n'},
        
        {0, 0, 0, 0}
    };

    while (opt != -1) {
        opt = getopt_long (argc, argv, "a:b:c:d:e:f:g:h:i:j:k:l:m:n", long_options, &opt_index);
        // cout << "opt is " << opt << "\n";
        switch(opt) {
            case 0:
            if (long_options[opt_index].flag != 0)                  break;
            case 'a':   // outfile_name             
                outfile_name = (string) optarg;                     break;
            case 'b':   // iyr              
                iyr = atoi(optarg);                                 break;
            case 'c':   // idoy             
                idoy= atoi(optarg);                                 break;
            case 'd':   // isec             
                isec= atoi(optarg);                                 break;
            case 'e':               
                flash_I0 = strtod(optarg, NULL);                    break;
            case 'f':   // lower rayfile                
                f1 = strtod(optarg, NULL);                          break;
            case 'g':   // upper rayfile                
                f2 = strtod(optarg, NULL);                          break; 
            case 'h':               
                ray_inp_dir = (string) optarg;                      break;
            case 'i':  // tmax              
                time_max = strtod(optarg, NULL);                    break;
            case 'j':  // max distance  
                max_ground_distance = strtod(optarg, NULL);         break;
            case 'k':  // num_times
                num_times = atoi(optarg);                           break;
            case 'l':  // frequency step size 
                num_freqs = strtod(optarg, NULL);                   break;
            case 'm':  // input latitude (geomagnetic)
                flash_pos[1] = strtod(optarg, NULL);                break;
            case 'n':  // input longitude (geomagnetic)
                flash_pos[2] = strtod(optarg, NULL);                break;
            case '?':
                 printf("\nUnknown option: %s\n",opt);  break;
        }
    }


    vector < Vector3d > prev_pos(num_freqs);
    vector < Vector3d > cur_pos(num_freqs);

    double time_step = (1.0*((1.0*time_max)/num_times));
    // int num_freqs = fmax(1,ceil((f2 - f1)/freq_step_size));

    // freq_step_size = fmin(freq_step_size, fabs(f2 - f1));
    double freq_step_size = fabs(f2 - f1)/(1.0*num_freqs);
    // // Storage space for the first sweep
    // vector<Vector3d> interp_points[num_freqs];
    // vector<double>   interp_data[num_freqs];

    double in_lat, in_lon, avg_distance_from_flash;
    map <int, rayF> raylist;
    rayF cur_rays[8];
    double latmin, latmax, lonmin, lonmax;
    double wmin, wmax, tmax;
    double dlat, dlon, dw;
    double lat0, lon0, rad0;
    double tmp_coords[3], flash_pos_sm[3];
    time_t run_tstart, run_tend;
    // int num_freqs_fine;
    double frame_area, initial_area;

    // Vector3d cell_pos;
    double path_length;

    Vector3d corners[8];
    double damping_at_corners[8];
    Vector3d centerpoint;
    double damping_avg;
    double bounding_sphere_radius;
    double kk;
    // Array of pointers to single-timestep frames
    rayT cur_frames[8];
    rayT prev_frames[8];


    int itime_in[2];
    int yearday = iyr*1000 + idoy;
    itime_in[0] = yearday;
    itime_in[1] = isec*1e3;


    // Talk about yaself
    cout << "---- power density calculation ---- " << endl;
    cout << "ray dir:\t" << ray_inp_dir << endl;
    cout << "outfile:\t" << outfile_name << endl;
    cout << "yearday:\t" << itime_in[0] << endl;
    cout << "msecs:\t" << itime_in[1] << endl;
    cout << "lat:\t" << flash_pos[1] << endl;
    cout << "lon:\t" << flash_pos[2] << endl;
    cout << "freqs:\t" << f1 << " to " << f2 << endl;
    cout << "step size:\t" << freq_step_size << " Hz" << endl;


// --------------- Get flash input coordinates -----------------------    
    lat0 = D2R*flash_pos[1];
    lon0 = D2R*flash_pos[2];
    rad0 = flash_pos[0];

    // Get flash position in SM coords:
    pol_to_cart_d_(&lat0, &lon0, &rad0, tmp_coords);
    mag_to_sm_d_(itime_in, tmp_coords, flash_pos_sm);


    ostringstream tmp;
    tmp << ray_inp_dir << "/f_" << f1;
    vector <vector<double> > available_rays;
    vector <vector<double> > adjacent_rays;
    vector <vector<double> > adjacent_rays_within_range;



    get_available_rays(tmp.str(), &available_rays);
    cout << "found " << available_rays.size() << " available rays\n";

    adjacent_rays = find_adjacent_rays(available_rays);
    cout << "found " << adjacent_rays.size() << " sets of adjacent rays\n";


    int num_rays_within_range = 0;

// -------------- Count how many adjacent ray sets we have within range: -----------------
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
        if ( (avg_distance_from_flash <= max_ground_distance) ) {
            adjacent_rays_within_range.push_back(*adj_row);
        }
    }

    int num_power_vects = adjacent_rays_within_range.size();
    cout << "num rays within range: " << num_power_vects << endl;


    // Allocate storage space for the power vectors:
    vector<Vector3d> interp_points[num_freqs][num_power_vects];
    vector<double>   energy_density[num_freqs][num_power_vects];
    vector<double>   geometric_spreading[num_freqs][num_power_vects];
    vector<double>   damping[num_freqs][num_power_vects];

// Loop back through the set of adjacent rays and calculate total input power: ----------

    // cout << "Integrating total input energy above iono:" << endl;
    // double total_integrated_power = total_input_power(flash_pos_sm, flash_I0, 
    //                                     15, 60, 70, 80, 200*2*PI, 30000*2*PI, itime_in);
    // cout << "total input energy: " << total_integrated_power << endl;


    // -------------- Iterate through adjacent ray sets: -----------------
    int set_index;
    for (vector < vector<double> >::iterator adj_row=adjacent_rays_within_range.begin(); adj_row!=adjacent_rays_within_range.end(); ++adj_row) {

        set_index = distance(adjacent_rays_within_range.begin(), adj_row);
        cout << "index: " << set_index << endl;
        avg_distance_from_flash = 0;

        for (int i=0; i<4; ++i) {
            in_lat = (*adj_row)[2*i  ];
            in_lon = (*adj_row)[2*i+1];
            avg_distance_from_flash += haversine_distance(in_lat, in_lon, flash_pos[1], flash_pos[2]);
        }
        avg_distance_from_flash /= 4000.0;  // Average dist of the 4 corner rays, in km
                                            // (Assuming both frequencies have same ray spacing) 
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

            cout << cur_rays[i].w/2/PI << " " << cur_rays[i].in_lat << " " << cur_rays[i].in_lon << " " << cur_rays[i].time.back() << endl;
        }
        
        tmax = min(tmax, time_max);

        // starting separation in lat, lon directions (meters)
        dlat = D2R*(R_E + H_IONO)*(latmax - latmin);
        dlon = D2R*(R_E + H_IONO)*(lonmax - lonmin)*cos(D2R*(latmax + latmin)/2.);
        dw   = wmax - wmin;

        cout << "\n----- current rays: -----\n";
        cout << "lon: " << lonmax << ", " << lonmin << " deg\n";
        cout << "lat: " << latmax << ", " << latmin << " deg\n";
        cout << "f: " << wmax/(2*PI) << ", " << wmin/(2*PI) << " Hz\n"; 
        cout << "Avg distance from flash: " << avg_distance_from_flash << " km\n";

        // Get total input power within each frequency band:
        double inp_pwr[num_freqs];
        for (int nf = 0; nf<num_freqs; ++nf) {
            double wlo = wmin + nf*dw/num_freqs;
            double whi = wmin + (nf + 1)*dw/num_freqs;
            inp_pwr[nf] = total_input_power(flash_pos_sm, flash_I0, 
                                        latmin, latmax, lonmin, lonmax, wlo, whi, itime_in);
            cout << "input power between " << wlo/2./PI << " and " << whi/2./PI << " : " << inp_pwr[nf] << " Joules\n";
        }

        // --------------------- Interpolate + sum over output grid ----------------
        //                                ( The main event)            
        // -------------------------------------------------------------------------
        time(&run_tstart);


        // Initial area (only needed if you're calculating relative spreading instead of absolute)
        initial_area = polygon_frame_area(corners);

        // Step forward in time:
        for (double tt=0; tt < tmax; tt+=time_step) {
            // interpolate current frames:
            for (int zz=0; zz<8; zz++) { interp_rayF(&cur_rays[zz], &(cur_frames[zz]), tt);}

            // for (double kk=0; kk < 1; kk += 1./num_freqs_fine) {
            for (int k_ind = 0; k_ind < num_freqs; ++k_ind) {

                kk = (1.0*k_ind)/(1.0*num_freqs);
                // cout << tt << " " << k_ind << " " << num_freqs << endl;


                // interpolate corners over frequency axis:
                for (int i=0; i<4; ++i) {
                    corners[i]   = cur_frames[i].pos*kk + cur_frames[i+4].pos*(1-kk);
                    // corners[i+4] = prev_frames[i].pos*kk + prev_frames[i+4].pos*(1-kk);
                    damping_at_corners[i] = cur_frames[i].damping*kk + cur_frames[i+4].damping*(1-kk);
                    // damping_at_corners[i+4] = prev_frames[i].damping*kk + prev_frames[i+4].damping*(1-kk);
                }
                
                // Average all four corners to get values at the center of the plane:
                cur_pos[k_ind][0] = 0; cur_pos[k_ind][1] = 0; cur_pos[k_ind][2] = 0;
                damping_avg = 0;
                for (int i=0; i<4; ++i) {
                    cur_pos[k_ind] += corners[i];
                    damping_avg += damping_at_corners[i];
                }
                cur_pos[k_ind] /=4.0;
                damping_avg /=4.0;


                // Get frame area for geometric factor:
                frame_area = polygon_frame_area(corners);

                // Get path length for energy density:
                if (tt==0) { 
                    path_length = 1; 
                } else {
                    path_length = (cur_pos[k_ind] - prev_pos[k_ind]).norm()*R_E;
                }
                
                interp_points[k_ind][set_index].push_back(cur_pos[k_ind]);

                //                            (    Joules /  m^2)      ( unitless ) / Meters  ~ Joules/m^3
                // double cur_energy_density = (inp_pwr[k_ind]/frame_area)*damping_avg/path_length;
                double cur_energy_density = (inp_pwr[k_ind])*damping_avg;   // <- total energy:  (so you can divide by the total volume in Python)

                energy_density[k_ind][set_index].push_back(cur_energy_density);

                geometric_spreading[k_ind][set_index].push_back(initial_area/frame_area);
                
                damping[k_ind][set_index].push_back(damping_avg);

                prev_pos[k_ind] = cur_pos[k_ind];
            }

            // // Step forward one frame:
            // for (int zz=0; zz<8; zz++) { prev_frames[zz] = cur_frames[zz]; }        
        }

        cout << "lengths of each  vector: ";
        for (int ii=0; ii<num_freqs; ii++) {
            cout << energy_density[ii][set_index].size() << " ";
        }
        cout << endl;

        }

        cout << "requested timesteps: " << num_times;
        // cout << " Actual timesteps: " << min_num_times << endl;
        FILE* outfile; 
        outfile = fopen(outfile_name.c_str(),"wb");
        if (outfile==NULL) {
            cout << "failed to open output file ;~;" << endl;
        } else {

            fprintf(outfile, "%d\t%d\t",num_freqs, num_power_vects);

            // Print lengths of each vector at the top of the file
            for (int k_ind = 0; k_ind < num_freqs; ++k_ind) {
                for (int i_ind = 0; i_ind < num_power_vects; i_ind++) {
                    fprintf(outfile, "%i ",energy_density[k_ind][i_ind].size());
                }
            }

            for (int k_ind = 0; k_ind < num_freqs; ++k_ind) {
                for (int i_ind = 0; i_ind < num_power_vects; i_ind++) {
                    for (int t_ind = 0; t_ind < energy_density[k_ind][i_ind].size(); ++t_ind) {
                        fprintf(outfile,"%g\t%g\t%g\t%g\t",
                            interp_points[k_ind][i_ind][t_ind][0],
                            interp_points[k_ind][i_ind][t_ind][1],
                            interp_points[k_ind][i_ind][t_ind][2],
                            energy_density[k_ind][i_ind][t_ind]);

                        fprintf(outfile,"%g\t", geometric_spreading[k_ind][i_ind][t_ind]);
                        fprintf(outfile,"%g\t", damping[k_ind][i_ind][t_ind]);
                    }
                }
            }
        }
        fclose(outfile);
}





    // ----------------- Scraps: Alglib stuff, nonlinear solver stuff... ------------








        // for (int k_ind = 0; k_ind < num_freqs; ++k_ind) {

        //     // // starting coords: (just print them for now)
        //     // for (int i_ind = 0; i_ind < num_power_vects; i_ind++){



        //     //     cout << interp_points[k_ind][i_ind][0].transpose() << endl;
        //     // }
        //     for (double tt=TIME_STEP; tt < tmax; tt+=TIME_STEP) {
        //         double cmin[3] = {1000, 1000, 1000};
        //         double cmax[3] = {-1000,-1000,-1000};
        //         // Find min and max indices:
        //         for (int i=0; i<num_power_vects; ++i) {
        //             for (int j=0; j<3; ++j) {
        //                 double tmp = interp_points[k_ind][i][tt-1].data()[j];
        //                 if (tmp > cmax[j]) {cmax[j] = tmp; }
        //                 if (tmp < cmin[j]) {cmin[j] = tmp; }
        //                 tmp = interp_points[k_ind][i][tt].data()[j];
        //                 if (tmp > cmax[j]) {cmax[j] = tmp; }
        //                 if (tmp < cmin[j]) {cmin[j] = tmp; }
        //             }
        //         }


        //         int xmin_ind = nearest(xaxis, NX, cmin[0], false);
        //         int xmax_ind = nearest(xaxis, NX, cmax[0], false);
        //         int ymin_ind = nearest(yaxis, NX, cmin[1], false);
        //         int ymax_ind = nearest(yaxis, NX, cmax[1], false);
        //         int zmin_ind = nearest(zaxis, NX, cmin[2], false);
        //         int zmax_ind = nearest(zaxis, NX, cmax[2], false); 

        //         cout << xmin_ind << " " << xmax_ind << " ";
        //         cout << ymin_ind << " " << ymax_ind << " ";
        //         cout << zmin_ind << " " << zmax_ind << endl;

        //         // And now you're back to the same problem
        //         // (but with correct center values).   

        //         // --> Decide how to interpolate over this set (2 x num_power_vects points).
        //         // This defines the wavefront at timestep tt.
        //         // Nearest neighbor / distance weighting? Area/volume weighting?
        //         // Fit a function to the wavefront positions and amplitudes? That's kinda cool...
                





            




        //     }
        // }


































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
    // FILE* outfile;
    // outfile = fopen(outfile_name.c_str(),"wb");
    // if (outfile==NULL) {
    //     cout << "failed to open output file ;~;" << endl;
    // } else {
    //     fwrite(out_grid, NX*NY*NZ*sizeof(double), 1, outfile);
    // }
    // fclose(outfile);

