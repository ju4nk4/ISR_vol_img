// 2024-07-01, Juan Carlos Araújo, ju4nk4@gmail.com

#include <getopt.h>
#include <stdio.h>
#include <sys/stat.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

//#include "interpolation_multi_parallel.hpp"
#include "omp.h"

//#include "src/ConvexHull.hpp"

#include "lac.hpp"
//#include "src/reader_hdf5.hpp"

#include "Ionosradar.hpp"

using namespace std;




// search for data:
// https://isr.sri.com/madrigal/cgi-bin/getMadfile.cgi?fileName=/opt/madrigal/experiments/2017/pfa/28jan17l/pfa170128.002

int main(int argc, char** argv) {
  auto start = std::chrono::system_clock::now();
  std::time_t end_time = std::chrono::system_clock::to_time_t(start);

  std::cout << "% Date: " << std::ctime(&end_time);
  printf(
      "%% Ionosradar: Creating a volumetric image from Incoherent Scatter "
      "Radar data.\n");

  unsigned int NUM_THREADS = 1, nx = 32, ny = 32, nz = 32,
               shepard_parameter = 3;  //, num_scalar_data = 5;

  bool estimate_derivatives = true, verbose = false;

  // int t0 = 0, tf = -1, test_case = 0;
  // bool specified_times = false;

  string directory = "",
         filename = "ras200113.004.hdf5";  //   /// ras151208.004 // "ras161121.002"; //

  const char* const short_options = "a:n:s:r:i:f:p:kvF:t:d:D:bV";
  // A string listing valid short options letters.

  double smothness = 0.5, ratio_R = 0.15;  //

  const struct option long_options[] = {
      // An array describing valid long options.
      {"alpha", required_argument, 0, 'a'},
      {"box_size", required_argument, 0, 'n'},
      {"smothness", required_argument, 0, 's'},
      {"ratio_r", required_argument, 0, 'r'},
      {"t0", required_argument, 0, 'i'},
      {"tf", required_argument, 0, 'f'},
      {"filename", required_argument, 0, 'F'},
      {"shepard_parameter", required_argument, 0, 'p'},
      {"estimate_derivatives", required_argument, 0, 'd'},
      {"data_directory", required_argument, 0, 'D'},
      {"num_threads", required_argument, 0, 't'},
      {"num_scalar_data", no_argument, 0, 'k'},
      {"beam_output", no_argument, 0, 'b'},
      {"plotting", no_argument, 0, 'V'},
      {NULL, 0, NULL, 0}};

  int next_option;
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, NULL);

    switch (next_option) {
      case 'n':  // -n or --box_size
        nx = atoi(optarg);
        ny = nx;
        nz = nx;
        break;
      case 's':  // -s or --smothness
        smothness = atof(optarg);
        break;
      case 'F':  // -F or --filename
        filename = (string)(optarg);
        // printf ("\n\n*********************** %s ************************\n",
        // filename.c_str() );
        break;
      case 'd':  // -d or --estimate_derivatives
        estimate_derivatives = atoi(optarg);
        break;
      case 'D':  // -D or --data_directory
        directory = (string)(optarg);
        break;
      case 't':    // --num_threads ot -t
        NUM_THREADS = atoi(optarg);
        break;
      case 'p':  // -p or --shepard_parameter
        shepard_parameter = atoi(optarg);
        break;
      case 'v':
        verbose = true;
      case '?':  // The user specified an invalid option.
        break;
      case -1:  // Done with options.
        break;
      default:  // Something else: unexpected.
        cout << "-- unexpected argument --" << endl;
        // abort ();
    }
  } while (next_option != -1);

  Ionosradar data;
  data.NUM_THREADS = NUM_THREADS;
  
  // Reading ISR data file
  // *************************************************************
  printf("%% Input:\n");
  if (directory != "") directory = directory + "/";

  const std::string path = directory + "data/" + filename;
  printf("%% \tDirectory: %s\n%% \tData file: %s\n", directory.c_str(),
         path.c_str());

  // Initializing Ionosradar
  // ***********************************************************
  if (smothness < 1e-10)
    smothness = 1e-10;  // Fixing wrong smoothness ... bug for s=0.0?
  data.set_interpolation_parameters(smothness, shepard_parameter, ratio_R,
                                    estimate_derivatives);  //

  printf("%% Loading data\n");
  data.load_radar_data_info(directory, filename);
  if (verbose) data.set_verbose();
  data.print_beam_info();
  data.compute_interpolation_mesh(nx, ny, nz);
  data.load_interpolate_Ion_in_time();  // ION, s: smoothness

  return 0;
}
