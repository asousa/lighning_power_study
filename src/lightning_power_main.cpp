#include <lightning_power.h>

int main(int argc, char *argv[]) 
{

    //                    alt, lat, lon (geomagnetic)
    double flash_pos[3] = {1, 50, 84};
    double flash_I0 = -10000;

    string ray_inp_dir = "/shared/users/asousa/WIPP/lightning_power_study/rays/globe_ngo";
    string outfile_name = "/shared/users/asousa/WIPP/lightning_power_study/rays/globe_ngo/test_dump.dat";


    double f1 = 200;
    double f2 = 230;
    int iyr = 2001;
    int idoy = 001;
    int isec = 0;

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
                for (double kk=0; kk < 1; kk += 1./num_freqs_fine) {


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

                    // centerpoint[0] = 0; centerpoint[1] = 0; centerpoint[2] = 0;
                    // // bounding sphere:
                    // bounding_sphere_radius = bounding_sphere(corners, centerpoint.data());
                    // // cout << " centerpoint: " <<  centerpoint.transpose() << " r: " << bounding_sphere_radius << endl;
                    // // bounding cube:
                    // int xmin_ind = nearest(xaxis, NX, centerpoint[0] - bounding_sphere_radius, false);
                    // int xmax_ind = nearest(xaxis, NX, centerpoint[0] + bounding_sphere_radius, false);
                    // int ymin_ind = nearest(yaxis, NY, centerpoint[1] - bounding_sphere_radius, false);
                    // int ymax_ind = nearest(yaxis, NY, centerpoint[1] + bounding_sphere_radius, false);
                    // int zmin_ind = nearest(zaxis, NZ, centerpoint[2] - bounding_sphere_radius, false);
                    // int zmax_ind = nearest(zaxis, NZ, centerpoint[2] + bounding_sphere_radius, false);

                    // Find min and max values to search thru:
                    double cmin[3] = {1000000, 1000000, 1000000};
                    double cmax[3] = {-1000000, -1000000, -1000000};
                    for (int zz=0; zz<8; ++zz) {
                        // cout << "corner " << zz << ":" << corners[zz].transpose() << endl;
                        for (int bb=0; bb<3; ++bb) {
                            if (corners[zz][bb] < cmin[bb]) {cmin[bb] = corners[zz][bb];}
                            if (corners[zz][bb] > cmax[bb]) {cmax[bb] = corners[zz][bb];}
                        }
                    }
                    // cout << "cmin: ";
                    // print_array(cmin,3);
                    // cout << "cmax: ";
                    // print_array(cmax,3);

                    int xmin_ind = int(floor((cmin[0] - XMIN)/GRID_STEP_SIZE));
                    int xmax_ind = int(floor((cmax[0] - cmin[0])/GRID_STEP_SIZE + xmin_ind));
                    int ymin_ind = int(floor((cmin[1] - YMIN)/GRID_STEP_SIZE));
                    int ymax_ind = int(floor((cmax[1] - cmin[1])/GRID_STEP_SIZE + ymin_ind));
                    int zmin_ind = int(floor((cmin[2] - ZMIN)/GRID_STEP_SIZE));
                    int zmax_ind = int(floor((cmax[2] - cmin[2])/GRID_STEP_SIZE + zmin_ind));

                    cout << xmin_ind << " " << xmax_ind << " ";
                    cout << ymin_ind << " " << ymax_ind << " ";
                    cout << zmin_ind << " " << zmax_ind << " " << endl;



                    // // ---------------- This block to just round to the nearest point
                    // for (int b = 0; b < 8; b++ ) {
                    //     x_ind = nearest(xaxis, NX, cur_frames[b].pos[0], false);
                    //     y_ind = nearest(yaxis, NX, cur_frames[b].pos[1], false);
                    //     z_ind = nearest(zaxis, NX, cur_frames[b].pos[2], false);

                    //     // x_ind = round((cur_frames[b].pos[0] - XMIN)/GRID_STEP_SIZE);
                    //     // y_ind = round((cur_frames[b].pos[1] - YMIN)/GRID_STEP_SIZE);
                    //     // z_ind = round((cur_frames[b].pos[2] - ZMIN)/GRID_STEP_SIZE);
                    //     out_grid[x_ind][y_ind][z_ind] +=  (inp_pwr/frame_area)*cur_frames[b].damping;

                    //     cout << x_ind << " " << y_ind << " " << z_ind << endl;
                    // }


                    if (xmin_ind >= 0 and xmax_ind < NX and
                        ymin_ind >= 0 and ymax_ind < NY and
                        zmin_ind >= 0 and zmax_ind < NZ) {
                        // Check for intercepts at each entry in the output space                    
                        for (int x_ind = xmin_ind; x_ind < xmax_ind; ++x_ind){
                            for (int y_ind = ymin_ind; y_ind < ymax_ind; ++y_ind){
                                for (int z_ind = zmin_ind; z_ind < zmax_ind; ++z_ind){

                                    // coordinates in center of cell:
                                    cell_pos[0] = XMIN + x_ind*GRID_STEP_SIZE + GRID_STEP_SIZE/2.;
                                    cell_pos[1] = YMIN + y_ind*GRID_STEP_SIZE + GRID_STEP_SIZE/2.;
                                    cell_pos[2] = ZMIN + z_ind*GRID_STEP_SIZE + GRID_STEP_SIZE/2.;


                                    // // Interpolate damping...
                                    // double wx = fabs((cmax[0] - cell_pos[0])/(cmax[0]-cmin[0]));
                                    // double wy = fabs((cmax[1] - cell_pos[1])/(cmax[1]-cmin[1]));
                                    // double wz = fabs((cmax[2] - cell_pos[2])/(cmax[2]-cmin[2]));

                                    // // cout << "weights: " << wx << " " << wy << " " << wz << endl;
                                    // double dd = 0;
                                    // dd += damping_at_corners[0]*(1-wx)*(1-wy)*(1-wz);
                                    // dd += damping_at_corners[1]*(wx)*(1-wy)*(1-wz);
                                    // dd += damping_at_corners[2]*(1-wx)*(wy)*(1-wz);
                                    // dd += damping_at_corners[3]*(wx)*(wy)*(1-wz);
                                    // dd += damping_at_corners[4]*(1-wx)*(1-wy)*(wz);
                                    // dd += damping_at_corners[5]*(wx)*(1-wy)*(wz);
                                    // dd += damping_at_corners[6]*(1-wx)*(wy)*(wz);
                                    // dd += damping_at_corners[7]*(wx)*(wy)*(wz);

                                    // out_grid[x_ind][y_ind][z_ind] += (inp_pwr/frame_area)*dd*FREQ_STEP_SIZE;
                                







                                    // out_grid[x_ind][y_ind][z_ind] += (inp_pwr/frame_area)*
                                    //                             interp_damping_inv_dist(corners, damping_at_corners, cell_pos) *
                                    //                             FREQ_STEP_SIZE;
                                

                                    find_crossing(corners, damping_at_corners, cell_pos, &ii, &jj);

                                    if (ii >=0 and ii <=1 and jj >= 0 and jj <= 1) {
                                        // cout << ii << " " << jj << endl;
                                        out_grid[x_ind][y_ind][z_ind] += (inp_pwr/frame_area);
                                    }            

                                }
                            }
                        }
                    }





                    // // Find min and max values to search thru:
                    // double cmin[3] = {1000};
                    // double cmax[3] = {-1000};
                    // for (int zz=0; zz<4; ++zz) {
                    //     for (int bb=0; bb<3; ++bb) {
                    //         if (cur_frames[zz].pos[bb]*kk + cur_frames[zz+4].pos[bb]*(1-kk)  < cmin[bb]) {
                    //             cmin[bb] = cur_frames[zz].pos[bb]*kk + cur_frames[zz].pos[bb]*(1 -kk);
                    //         }
                    //         if (cur_frames[zz].pos[bb]*kk + cur_frames[zz+4].pos[bb]*(1-kk) > cmax[bb]) {
                    //             cmax[bb] = cur_frames[zz].pos[bb]*kk + cur_frames[zz+4].pos[bb]*(1-kk);
                    //         }
                    //     }
                    // }

                    // int xmin_ind = nearest(xaxis, NX, cmin[0], false);
                    // int xmax_ind = nearest(xaxis, NX, cmax[0], false);
                    // int ymin_ind = nearest(yaxis, NY, cmin[1], false);
                    // int ymax_ind = nearest(yaxis, NY, cmax[1], false);
                    // int zmin_ind = nearest(zaxis, NZ, cmin[2], false);
                    // int zmax_ind = nearest(zaxis, NZ, cmax[2], false);


                    // // ---------------- This block to just round to the nearest point
                    // for (int b = 0; b < 8; b++ ) {
                    //     x_ind = nearest(xaxis, NX, cur_frames[b].pos[0], false);
                    //     y_ind = nearest(yaxis, NX, cur_frames[b].pos[1], false);
                    //     z_ind = nearest(zaxis, NX, cur_frames[b].pos[2], false);

                    //     // x_ind = round((cur_frames[b].pos[0] - XMIN)/GRID_STEP_SIZE);
                    //     // y_ind = round((cur_frames[b].pos[1] - YMIN)/GRID_STEP_SIZE);
                    //     // z_ind = round((cur_frames[b].pos[2] - ZMIN)/GRID_STEP_SIZE);
                    //     out_grid[x_ind][y_ind][z_ind] +=  (inp_pwr/frame_area)*cur_frames[b].damping;

                    //     cout << x_ind << " " << y_ind << " " << z_ind << endl;
                    // }

                    // cout << "kk: " << kk << " tt: " << tt << endl;


















                    // if (xmin_ind >= 0 and xmax_ind < NX and
                    //     ymin_ind >= 0 and ymax_ind < NY and
                    //     zmin_ind >= 0 and zmax_ind < NZ) {
                    //     // Check for intercepts at each entry in the output space                    
                    //     for (int x_ind = xmin_ind; x_ind < xmax_ind; ++x_ind){
                    //         for (int y_ind = ymin_ind; y_ind < ymax_ind; ++y_ind){
                    //             for (int z_ind = zmin_ind; z_ind < zmax_ind; ++z_ind){

                    //                 // coordinates in center of cell:
                    //                 cell_pos[0] = XMIN + x_ind*GRID_STEP_SIZE;
                    //                 cell_pos[1] = YMIN + y_ind*GRID_STEP_SIZE;
                    //                 cell_pos[2] = ZMIN + z_ind*GRID_STEP_SIZE;

                    //                 out_grid[x_ind][y_ind][z_ind] += (inp_pwr/frame_area)*
                    //                                                 interp_damping_inv_dist(cur_frames, prev_frames, kk, cell_pos) *
                    //                                                 FREQ_STEP_SIZE;
                    //                 // // cout << out_grid[x_ind][y_ind][z_ind] << endl;
                    //                 // find_crossing(cur_frames, prev_frames, cell_pos, kk, &ii, &jj);            

                    //                 // // check if crossing is within bounds:
                    //                 // if ( (ii >= 0 ) && (ii <= 1) && (jj >=0) && (jj <=1)) {
                    //                 //     out_grid[x_ind][y_ind][z_ind] += interp_damping(cur_frames, ii, jj, kk);
                    //                 //     // cout << "tt: " << tt << " ii: " << ii << " jj: " << jj << " kk: " << kk << endl;
                    //                 // }
                    //             }
                    //         }
                    //     }
                    // }
                }
            // Step forward one frame:
            for (int zz=0; zz<8; zz++) { prev_frames[zz] = cur_frames[zz]; }        }

            }
        }


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