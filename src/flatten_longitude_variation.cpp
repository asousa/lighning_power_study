#include <lightning_power.h>
// #include "stdafx.h"
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include "interpolation.h"


int main(int argc, char *argv[]) 
{



    string inp_dir = "/shared/users/asousa/WIPP/rays/2d/nightside/gcpm_kp0";
    string out_dir = "/shared/users/asousa/WIPP/rays/2d/nightside/flat/gcpm_kp0";


    int iyr = 2010;
    int idoy = 001;
    int isec = 0;
    

    int opt;
    int opt_index;
    static struct option long_options[] =
    {
        {"inp_dir",       required_argument,    0, 'a'},
        {"out_dir",       required_argument,    0, 'b'},
        {"iyr",           required_argument,    0, 'c'},
        {"idoy",          required_argument,    0, 'd'},
        {"isec",          required_argument,    0, 'e'},
        {0, 0, 0, 0}
    };

    while (opt != -1) {
        opt = getopt_long (argc, argv, "a:b:c:d:e:", long_options, &opt_index);
        // cout << "opt is " << opt << "\n";
        switch(opt) {
            case 0:
            if (long_options[opt_index].flag != 0)                  break;
            case 'a':   // input directory             
                inp_dir = (string) optarg;                          break;
            case 'b':   // output directory             
                out_dir = (string) optarg;                          break;
            case 'c':   // iyr              
                iyr = atoi(optarg);                                 break;
            case 'd':   // idoy             
                idoy= atoi(optarg);                                 break;
            case 'e':   // isec             
                isec= atoi(optarg);                                 break;
            case '?':
                 printf("\nUnknown option: %s\n",opt);  break;
        }
    }

    map <int, rayF> raylist;


    int itime_in[2];
    int yearday = iyr*1000 + idoy;
    itime_in[0] = yearday;
    itime_in[1] = isec*1e3;

    vector <vector<double> > available_rays;
    ostringstream cur_file;
    struct stat sb;



    // Talk about yaself
    cout << "---- Longitude Variation Flattening Machine ---- " << endl;
    cout << "input dir:\t" << inp_dir << endl;
    cout << "output dir:\t" << out_dir << endl;
    cout << "yearday:\t" << itime_in[0] << endl;
    cout << "msecs:\t" << itime_in[1] << endl;

// --------------- Get flash input coordinates -----------------------    

    get_available_rays(inp_dir, &available_rays);
    cout << "found " << available_rays.size() << " available rays\n";


    for (vector < vector<double> >::iterator row=available_rays.begin(); row!=available_rays.end(); ++row) {
        ostringstream cur_file, cur_out_dir, out_file;
        rayF ray;
        print_vector(*row);

        int f   = (*row)[0];
        int lat = (*row)[1];
        int lon = (*row)[2];

        cur_file << inp_dir   << "/f_" << f   << "/lon_"    << lon
                 << "/ray_" << f  << "_"   << lat << "_" << lon << ".ray";

        cur_out_dir << out_dir   << "/f_" << f   << "/lon_"    << lon;
                 
        out_file << cur_out_dir.str() << "/ray_" << f  << "_"   << lat << "_" << lon << ".ray";

        if (stat(cur_out_dir.str().c_str(), &sb) == 0 &(sb.st_mode)){
        } else {
            // Out directory doesn't exist, so make it:
            cout << "making output directory " << cur_out_dir.str() << endl;
            ostringstream cmd;
            cmd << "mkdir -p " << cur_out_dir.str();
            system(cmd.str().c_str());
        }

        cout << cur_file.str() << endl;
        cout << out_file.str() << endl;
        raylist = read_rayfile(cur_file.str());

        // My rays are one ray per file, so I'm not gonna iterate over the whole set.
        ray = raylist.at(1); // First ray should have ray number "1"
        for (int i=0; i < ray.pos.size(); ++i) {
            double tmp_mag[3], tmp_sm[3];

            sm_to_mag_d_(itime_in, ray.pos[i].data(),tmp_mag);
            cardeg(tmp_mag);
            tmp_mag[2] = lon;
            degcar(tmp_mag);
            mag_to_sm_d_(itime_in, tmp_mag, ray.pos[i].data());
        }

        raylist[1] = ray;
        write_rayfile(out_file.str(), raylist);

    }
  
}





















