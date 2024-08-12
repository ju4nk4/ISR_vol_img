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


#include "lac.hpp"
#include "reader_hdf5.hpp"
#include "ConvexHull.hpp"


/*

H5::H5File
H5::Group
H5::DataSpace
H5::DataSet
H5::PredType

*/



class Ionosradar {
 private:
  std::string ROOT_GROUP, PAR_1D, PAR_2D,
      ANG_AZM_DATA,  // azimuth
      ANG_ELM_DATA,  // elevation
      RANGE_DATA, TIME_DATA, BEAM_DATA_NEL, BEAM_DATA_TI, BEAM_DATA_TE,
      BEAM_DATA_POp, BEAM_DATA_VO, BEAM_DATA_UNEL;

  // interpolation parameters
  double smoothness, ratio_R;
  int shepard_parameter;

  // Fourier coefficients for Convex-Hull
  std::vector<double> xA, xB, yA, yB;

 public:
  std::unique_ptr<H5::H5File> _file = nullptr;
  
  const double Re = 6378.137; 
  unsigned int NUM_THREADS = 1;
  
  // Instrument data
  double ISR_lat, ISR_lon, ISR_alt, Re_alt, min_alt, max_alt, res_z = 1e-10;
  Metadata ISR_meta;
  std::string start_time, end_time;

  std::vector<H5::Group> beamGroups;
  H5::Group beam_closest_z_axis;
  std::vector<BeamInfo> beam_info;  // beam_data

  std::array<double, 6> data_bounding_box = {1e8, -1e8, 1e8, -1e8, 1e8, -1e8};
  bool estimate_derivatives = true, verbose = false;

  std::vector<glm::dvec3> unit_beam_dir_array;
  std::vector<std::array<double, 2>> min_max_range_array;
  glm::dvec3 avg_dir, unit_beam_test_array;
  glm::dvec3 M_mom;

  // Reference vectors in K cell
  glm::dvec3 K_ux, K_uy, K_uz;

  glm::dmat3 R_glob2loc;

  // Prepare for interpolation
  std::vector<std::array<double, 3>> P_data, P_max_data, P_max_polar;  //, P;

  // interpolation grid data
  unsigned int num_scalar_data, nx, ny, nz, ni,
      time_measurements;  // time_measurements = numOfDimSamples.rows
  int t0 = 0, tf = -1;
  double* xi;  // interpolating points

  std::vector<std::array<double, 3>> pi;

  double R, RM;  // interpolation limits

  std::fstream file;
  std::string data_name, data_path, output_folder, vtk_filename;
  
  // -------------------------------------------------------------------------------------
  void print_metadata();
  void set_verbose();
  void load_experiment_parameters(const std::string& path);
  void set_interpolation_parameters(double e_smoothness,
                                    int e_shepard_parameter, double e_ratio_R,
                                    bool e_estimate_derivatives);
  // -------------------------------------------------------------------------------------
  void load_radar_data_info(const std::string& directory,
                            const std::string& filename);
  // -------------------------------------------------------------------------------------
  void print_beam_info();
  void compute_interpolation_mesh(int e_nx, int e_ny, int e_nz);
  // -------------------------------------------------------------------------------------
  void load_ISR_data_in_time(std::vector<std::array<double, 8>>& radar_data,
                             int time_sample);
  // -------------------------------------------------------------------------------------
  void load_interpolate_Ion_in_time();
  // -------------------------------------------------------------------------------------
  void set_bounding_box(std::array<double, 6>& in_bounding_box);
  
  double min_angle_sq_fit_data(std::vector<std::array<double, 2>>& X);
  
  void compute_fourier_convex_hull(std::vector<std::array<double, 2>>& coords,
                                   std::vector<double>& angles, double L_scale);

  void write_interpolated_time_vals_vtk(std::string data_name, unsigned int time,
                                        double* zi);

  void output_vtk();

  void output_radar_data_vtk(
      std::string name,
      std::vector<std::array<double, 3>>& P);
      
  void write_interpolated_time_selected_vals_vtk(
      std::string data_name, unsigned int n_data, unsigned int time,
      std::vector<std::array<double, 5>>& data_interp);
  // -------------------------------------------------------------------------------------
  void generateInterpolationTestGrid(int test_case);

  void generateInterpolationGrid(int e_nx, int e_ny, int e_nz);
  // -------------------------------------------------------------------------------------
  // Shepard interpolation for Ion data
  void sparseInterpolationAlgorithm_large(
      int timestep, std::vector<std::array<double, 8>>& data_array,
      std::string fname);
      
  double* shepard_interp_der_2d(int p, unsigned int nd, double* xd, double* fd,
                                unsigned int ni, double* xi,
                                bool estimate_derivatives, double R, double RM,
                                double c_nu);
  
  
 private:
  int _maxTime = 1;

  double  fourier_radius(double t);
  
  double* fourier_radius_xy(double t);

  glm::dvec3 fourier_xy(double t);

  double top_lead_interpolation(double x, double y, std::vector<glm::dvec3>& V,
                                unsigned int nd_e, double* xd_e,
                                double* fd_e);

   std::array<double, 3> transfinite_interpolation_Fourier_CH_Poly_top(
      double x, double y, double z, std::vector<glm::dvec3>& V,
      unsigned int n_poly, std::vector<double>& top_poly_points,
      std::vector<double>& top_poly_coeff, double R_rm, double R_RM,
      double LM);

  
  void scaled3D_TI_Fourier_CH_Poly_top_grid (
      unsigned int nx, unsigned int ny, unsigned int nz,
      std::vector<glm::dvec3>& V,
      // unsigned int nd_e, double *xd_e, double *fd_e,
      unsigned int n_poly, std::vector<double>& top_poly_points,
      std::vector<double>& top_poly_coeff,
      std::vector<std::array<double, 3>>& P, double R_rm, double R_RM,
      double LM);
  
  // -------------------------------------------------------------------------------------
  
  void collect_BeamData (
      H5::Group& beamGroup, std::vector<std::array<double, 3>>& points_array,
      std::vector<std::array<double, 3>>& points_max_array,
      std::vector<std::array<double, 3>>& points_max_array_polar,
      std::vector<glm::dvec3>& unit_beam_dir_array,
      std::vector<std::array<double, 2>>& min_max_range_array,
      std::array<double, 6>& beam_bounding_box);
  
  void collect_BeamData_Large(H5::Group& beamGroup,
                              std::vector<std::array<double, 8>>& data_array,
                              unsigned int time_step);
  
  
  // -----------------------------------------------------------------------------------
};  // Ionosradar






