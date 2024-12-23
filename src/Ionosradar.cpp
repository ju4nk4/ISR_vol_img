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

#include "interpolation_multi_parallel.hpp"
//#include "omp.h"

//#include "src/ConvexHull.hpp"

//#include "src/lac.hpp"
//#include "src/reader_hdf5.hpp"



/*

H5::H5File
H5::Group
H5::DataSpace
H5::DataSet
H5::PredType

*/

// Threshold NAN in ne
float TH_NAN_NE = 0.5f;

#include "Ionosradar.hpp"

using namespace std;


void Ionosradar::print_metadata() { ISR_meta.print(); }
void Ionosradar::set_verbose() { verbose = true; }
void Ionosradar::load_experiment_parameters(const std::string& path) {
    ROOT_GROUP = "Data/Array Layout";
    PAR_1D = "1D Parameters";
    PAR_2D = "2D Parameters";
    ANG_AZM_DATA = "azm";  // azimuth
    ANG_ELM_DATA = "elm";  // elevation
    RANGE_DATA = "range";
    TIME_DATA = "timestamps";
    BEAM_DATA_NEL = "nel";
    BEAM_DATA_TI = "ti";
    BEAM_DATA_TE = "te";
    BEAM_DATA_POp = "po+";
    BEAM_DATA_VO = "vo";
    BEAM_DATA_UNEL = "popl";
    num_scalar_data = 5;

    read_metadata(path, ISR_meta);
    print_metadata();

    double e_lat = std::stod(ISR_meta.inst_lat),
           e_lon = std::stod(ISR_meta.inst_lon),
           e_alt = std::stod(ISR_meta.inst_alt), deg_to_rad = M_PI / 180.0;

    ISR_lat = e_lat * deg_to_rad;
    ISR_lon = e_lon * deg_to_rad;
    ISR_alt = e_alt;

    printf("%% \tEntered: load_experiment_parameters. Opening: %s\n",
           path.c_str());

    Re_alt = Re + ISR_alt;  // Earth radius at sea level + instrument altitude

    printf("%% \tReading Metadata: ISR (lat, lon, alt), ti, tf. \n");

    // transforming the data
    double x_angle = 0.5 * M_PI - ISR_lat;
    double z_angle = 0.5 * M_PI + ISR_lon;

    glm::dmat3 iRotX = MatRotateAxis(-x_angle, glm::dvec3(1.0, 0.0, 0.0));
    glm::dmat3 iRotZ = MatRotateAxis(-z_angle, glm::dvec3(0.0, 0.0, 1.0));
    R_glob2loc = matMul(iRotX, iRotZ);  // iRotX * iRotZ;

    printf("%% \tComputing global data: coord transformations. \n");
}

void Ionosradar::set_interpolation_parameters (double e_smoothness,
                                    int e_shepard_parameter, double e_ratio_R,
                                    bool e_estimate_derivatives) {
    smoothness = e_smoothness;
    shepard_parameter = e_shepard_parameter;
    ratio_R = e_ratio_R;
    estimate_derivatives = e_estimate_derivatives;
}
// -------------------------------------------------------------------------------
  void Ionosradar::load_radar_data_info(const std::string& directory,
                            const std::string& filename) {
    data_name = filename;
    std::string path = directory + "data/" + filename;
    data_path = path;
    string spath = filename;

    load_experiment_parameters(path);

    printf("%% \tIonosradar: load \n");
    fflush(stdout);

    const H5std_string h5filename(path);
    try {
      _file = std::make_unique<H5::H5File>(path.c_str(), H5F_ACC_RDONLY);
    } catch (const H5::FileIException& err) {
      std::cerr << "Ionosradar: Error reading HDF5 data: {} " << std::endl;
    }

    std::replace(spath.begin(), spath.end(), '/', '_');
    std::replace(spath.begin(), spath.end(), '.', '_');

    output_folder = directory + "output";
    mkdir(output_folder.c_str(), 0777);

    output_folder = directory + "output/" + spath;
    mkdir(output_folder.c_str(), 0777);

    H5::Group data = _file->openGroup(ROOT_GROUP);
    beamGroups = getAllGroups(data);

    std::array<double, 6> data_bounding_box = {1e8, -1e8, 1e8, -1e8, 1e8, -1e8};

    for (H5::Group& beamGroup : beamGroups) {
      std::array<double, 6> beam_bounding_box;
      collect_BeamData(beamGroup, P_data, P_max_data, P_max_polar,
                       unit_beam_dir_array, min_max_range_array,
                       beam_bounding_box);

      data_bounding_box[0] = min(beam_bounding_box[0], data_bounding_box[0]);
      data_bounding_box[1] = max(beam_bounding_box[1], data_bounding_box[1]);

      data_bounding_box[2] = min(beam_bounding_box[2], data_bounding_box[2]);
      data_bounding_box[3] = max(beam_bounding_box[3], data_bounding_box[3]);

      data_bounding_box[4] = min(beam_bounding_box[4], data_bounding_box[4]);
      data_bounding_box[5] = max(beam_bounding_box[5], data_bounding_box[5]);
    }

    set_bounding_box(data_bounding_box);
    unsigned int acc = 0;
    for (unsigned int i = 0; i < beam_info.size(); i++) {
      beam_info[i].index_i = acc;
      beam_info[i].index_f = acc + beam_info[i].num_range - 1;

      acc = acc + beam_info[i].num_range;
    }

    printf("%% \tIonosradar: load_radar_data_info: \n\tNp = %lu; Nb = %lu;\n",
           P_data.size(), beam_info.size());
    fflush(stdout);

    t0 = 0;
    tf = time_measurements;
}
// ---------------------------------------------------------------------------------------
void Ionosradar::print_beam_info() {
    printf(
        "%%\tPrinting beam data info [beam_id, data ordering "
        "(index_i,index_f), num range, elm, azm ]:\n");
    printf("\tbeam_info = [\n");
    for (unsigned int i = 0; i < beam_info.size(); i++) {
      printf("\t\t%d, %3d, %3d, %d,% .3f,% .3f \n", beam_info[i].id,
             beam_info[i].index_i, beam_info[i].index_f, beam_info[i].num_range,
             beam_info[i].elm, beam_info[i].azm);
    }
    printf("\t];\n");
}
// ---------------------------------------------------------------------------------------
void Ionosradar::compute_interpolation_mesh(int e_nx, int e_ny, int e_nz) {
    nx = e_nx;
    ny = e_ny;
    nz = e_nz;

    printf("%% \tIonosradar: compute_interpolation_mesh \n");
    fflush(stdout);

    vtk_filename = "NEL";

    generateInterpolationGrid(e_nx, e_ny, e_nz);

    printf("%% \tIonosradar: compute_initial_phys_at_interp_nodes \n");
    fflush(stdout);
}
// ---------------------------------------------------------------------------------------
void Ionosradar::load_ISR_data_in_time(std::vector<std::array<double, 8>>& radar_data,
                             int time_sample) {
    for (H5::Group& beamGroup : beamGroups) {
      collect_BeamData_Large(beamGroup, radar_data, time_sample);
    }
  }
// ---------------------------------------------------------------------------------------
void Ionosradar::load_interpolate_Ion_in_time() {
    std::vector<std::array<double, 8>> radar_data, z_beam_data;

    for (int t = t0; t < tf; t++) {  // time samples

      collect_BeamData_Large(beam_closest_z_axis, z_beam_data, t);
      load_ISR_data_in_time(radar_data, t);
      float ratio_NAN_ne = ratioNANne(radar_data);
      if (ratio_NAN_ne < TH_NAN_NE) { 
        sparseInterpolationAlgorithm_large(t, radar_data, "test_conductivities");
      }
      else {
        produceNANResult(t);
      }

      radar_data.clear();
      z_beam_data.clear();
    }
  }
// ---------------------------------------------------------------------------------------
void Ionosradar::set_bounding_box(std::array<double, 6>& in_bounding_box) {
    printf("%% \tIonosradar: set_bounding_box (rectangular box from data)\n");
    printf("\tdata_bounding_box = [");
    for (unsigned int i = 0; i < 6; i++) {
      data_bounding_box[i] = in_bounding_box[i];
      printf("% .4f ", data_bounding_box[i]);
    }
    printf("\t];\n");
}
// -------------------------------------------------------------------------------------
  /*
  % Routine: projects/best_sq_outline/main.m
  % Try to find the best rotation of a square so that the radar beams (projected
  to a plane cut) % fit nicely inside. We want to parametrize the faces of the
  square to have a polar shape. % Good for AMISR data, but fails for arbitrary
  randomly generated data! % For the minimization, we use a flattened Gaussian
  function resembling a curved square. % Then, we rotate the wire (Convex Hull)
  until it falls by gravity (minimum). % The final rotation is the output.
  */

  double Ionosradar::min_angle_sq_fit_data(vector<std::array<double, 2>>& X) {
    unsigned int Np = X.size(), p = 6, N_sweep = 91;  // -45,0,45

    double r = 0.0, x, y, l, c, obj_f, min_f = 1.0e14, min_angle,
           c_rad = M_PI / 180, ang;

    for (unsigned int i = 0; i < Np; i++) {
      x = X[i][0];
      y = X[i][1];
      r += sqrt(x * x + y * y);
    }

    l = r / Np;  // mean r
    c = 0.5 * M_PI / l;

    for (unsigned int k = 0; k < N_sweep; k++) {
      ang = (-45.0 + k) * c_rad;

      double co = cos(ang), so = sin(ang), rp;

      obj_f = 0.0;

      for (unsigned int i = 0; i < Np; i++) {
        x = X[i][0] * co + X[i][1] * so;
        y = -X[i][0] * so + X[i][1] * co;

        rp = pow(c * x, p) + pow(c * y, p);
        obj_f += exp(-rp);
      }
      
      if (obj_f < min_f) {
        min_f = obj_f;
        min_angle = ang;
      }
    }

    printf(
        "%% \tIonosradar::min_angle_sq_fit_data: \n\tmin_ang = % f; min_f = % "
        "f;\n\n",
        (min_angle / c_rad), min_f);
    return (-1.0 * min_angle);
  }
// -----------------------------------------------------------------------------------
  // Compute Fourier representation of the Convex-Hull
  void Ionosradar::compute_fourier_convex_hull(vector<std::array<double, 2>>& coords,
                                   vector<double>& angles, double L_scale) {
    unsigned int Nf = 10, Np = angles.size() + 1;

    vector<double> Px(Np), Py(Np), o(Np);

    for (unsigned int i = 0; i < angles.size(); i++) {
      Px[i] = coords[i][0];
      Py[i] = coords[i][1];
      o[i] = angles[i];
    }

    Px[Np - 1] = Px[0];
    Py[Np - 1] = Py[0];
    o[Np - 1] = o[0] + 2.0 * M_PI;

    double xA0 = 0.0, yA0 = 0.0;

    xA.resize(Nf + 1);
    xB.resize(Nf + 1);
    yA.resize(Nf + 1);
    yB.resize(Nf + 1);

    for (unsigned int n = 0; n <= Nf; n++) {
      xA[n] = 0.0;
      yA[n] = 0.0;
      xB[n] = 0.0;
      yB[n] = 0.0;
    }

    for (unsigned int i = 0; i < Np - 1; i++) {  // % increased
      double a, b, xa, xb, xc, xm, ya, yb, yc, ym;

      a = o[i];
      b = o[i + 1];
      xa = Px[i];
      xb = Px[i + 1];
      ya = Py[i];
      yb = Py[i + 1];

      xm = (xb - xa) / (b - a);
      ym = (yb - ya) / (b - a);

      xc = xa - a * (xb - xa) / (b - a);
      yc = ya - a * (yb - ya) / (b - a);

      xA0 =
          xA0 + 0.25 * (b - a) * (xm * (a + b) + 2 * xc) / M_PI;  
      yA0 =
          yA0 + 0.25 * (b - a) * (ym * (a + b) + 2 * yc) / M_PI;  

      for (int n = 1; n <= (int)Nf; n++) {
        int n2 = n * n;
        double sna, cna, snb, cnb;

        sna = sin(n * a);
        cna = cos(n * a);
        snb = sin(n * b);
        cnb = cos(n * b);

        xA[n] = xA[n] + ((-n * (a * xm + xc) * sna - xm * cna +
                          n * (b * xm + xc) * snb + xm * cnb) /
                         n2) /
                            M_PI;
        yA[n] = yA[n] + ((-n * (a * ym + yc) * sna - ym * cna +
                          n * (b * ym + yc) * snb + ym * cnb) /
                         n2) /
                            M_PI;

        xB[n] = xB[n] + ((xm * (snb - sna) + n * (a * xm + xc) * cna -
                          n * (b * xm + xc) * cnb) /
                         n2) /
                            M_PI;
        yB[n] = yB[n] + ((ym * (snb - sna) + n * (a * ym + yc) * cna -
                          n * (b * ym + yc) * cnb) /
                         n2) /
                            M_PI;
      }
    }

    xA[0] = xA0;
    yA[0] = yA0;
    xB[0] = 0.0;
    yB[0] = 0.0;

    // Normalize: we use the rate
    double xN = L_scale, yN = L_scale;  
    for (unsigned int n = 0; n <= Nf; n++) {
      xA[n] = xA[n] / xN;
      yA[n] = yA[n] / yN;

      xB[n] = xB[n] / xN;
      yB[n] = yB[n] / yN;
    }
  }
// -----------------------------------------------------------------------------------
  void Ionosradar::write_interpolated_time_vals_vtk(string data_name, unsigned int time,
                                        double* zi) {
    std::ostringstream oss;
    oss << std::setw(5) << std::setfill('0') << time;

    string filename = output_folder + "/" + "N" + std::to_string(nx) + "_" +
                      vtk_filename + "_" + oss.str() + ".vtk";

    printf("%% output: %s, %u out of %u\n", filename.c_str(), time,
           time_measurements);

    fstream file;
    file.open(filename.c_str(), ios::out);

    unsigned int n = nx * ny * nz;  //, ni = pi.size();

    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    file << "POINTS " << n << " float\n";

    for (unsigned int i = 0; i < ni; i++) {
      file << (pi[i][0]) << " " << (pi[i][1]) << " " << (pi[i][2]) << "\n";
    }

    file << "\n";
    file << "POINT_DATA " << ni << "\n";
    file << "SCALARS " << "nel" << " float\n";
    file << "LOOKUP_TABLE default\n";

    for (unsigned int i = 0; i < ni; i++) {
      file << (zi[i]) << " ";
    }
    file << "\n";

    file.close();
  }

  void Ionosradar::output_vtk() { file.close(); }

// -----------------------------------------------------------------------------------
  void Ionosradar::output_radar_data_vtk(
      string name,
      std::vector<std::array<double, 3>>& P) {

    string filename = output_folder + "/" + name, file_vtk;

    printf("%% output: %s.vtk\n", filename.c_str());

    fstream file_data;
    file_vtk = filename + ".vtk";
    file_data.open(file_vtk.c_str(), ios::out);

    unsigned int n = P.size();

    file_data << "# vtk DataFile Version 3.0\n";
    file_data << "vtk output\n";
    file_data << "ASCII\n";
    file_data << "DATASET POLYDATA\n";
    file_data << "POINTS " << n << " float\n";

    for (unsigned int i = 0; i < n; i++) {
      file_data << (P[i][0]) << " " << (P[i][1]) << " " << (P[i][2]) << "\n";
    }

    file_data << "\n";
    file_data.close();
  }
// -----------------------------------------------------------------------------------
void Ionosradar::write_interpolated_time_selected_vals_vtk(
      string data_name, unsigned int n_data, unsigned int time,
      std::vector<std::array<double, 5>>& data_interp) {
    std::ostringstream oss;
    oss << std::setw(5) << std::setfill('0') << time;

    string filename = output_folder + "/" + "N" + std::to_string(nx) + "_" +
                      data_name + "_" + oss.str() + ".vtk";

    printf("output: %s, out of %u\n", filename.c_str(), time_measurements);

    fstream file;
    file.open(filename.c_str(), ios::out);

    unsigned int n = nx * ny * nz;

    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    file << "POINTS " << n << " float\n";

    for (unsigned int i = 0; i < ni; i++) {
      file << (float)(pi[i][0]) << " " << (float)(pi[i][1]) << " "
           << (float)(pi[i][2]) << "\n";
    }
    file << "\n";

    file << "POINT_DATA " << ni << "\n";

    // -------------------------------------------------------------------------------
    file << "SCALARS " << "nel" << " float\n";
    file << "LOOKUP_TABLE default\n";

    for (unsigned int i = 0; i < ni; i++) {
      file << (float)(data_interp[i][0]) << " ";  
    }
    file << "\n";

    // -------------------------------------------------------------------------------
    file << "SCALARS " << "nOp" << " float\n";
    file << "LOOKUP_TABLE default\n";

    for (unsigned int i = 0; i < ni; i++) {
      file << (float)(data_interp[i][1]) << " ";  //
    }
    file << "\n";

    // -------------------------------------------------------------------------------
    file << "SCALARS " << "Ti" << " float\n";
    file << "LOOKUP_TABLE default\n";

    for (unsigned int i = 0; i < ni; i++) {
      file << (float)(data_interp[i][3]) << " ";  
    }
    file << "\n";

    // -------------------------------------------------------------------------------
    file << "SCALARS " << "Te" << " float\n";
    file << "LOOKUP_TABLE default\n";

    for (unsigned int i = 0; i < ni; i++) {
      file << (float)(data_interp[i][4]) << " ";  
    }
    file << "\n";

    // -------------------------------------------------------------------------------
    file << "SCALARS " << "Vlos" << " float\n";
    file << "LOOKUP_TABLE default\n";

    for (unsigned int i = 0; i < ni; i++) {
      file << (float)(data_interp[i][2]) << " ";  
    }
    file << "\n";

    file << "\n";
    file.close();
  }
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
  void Ionosradar::generateInterpolationTestGrid(int test_case) {
    if (test_case == 0) {
      nx = 0;
      ny = 0;
      nz = 1000;

      string spoints = "scattered_points";  
      output_radar_data_vtk(
          spoints,
          P_data);  

      ni = nz;
      xi = new double[3 * ni];  // Global in the class

      double r0 = min_alt, r1 = max_alt, dr = (r1 - r0) / (ni - 1);

      glm::dvec3 pos = glm::dvec3(0, 0, 0);

      // Global RM
      RM = (r1 - r0);

      for (unsigned int i = 0; i < ni; i++) {
        pos = (r0 + i * dr) * unit_beam_test_array;

        xi[0 + 3 * i] = pos.x;
        xi[1 + 3 * i] = pos.y;
        xi[2 + 3 * i] = pos.z;

        pi.push_back({xi[0 + 3 * i], xi[1 + 3 * i], xi[2 + 3 * i]});
      }
    }  
  }
// -----------------------------------------------------------------------------------
// when we have the first time sample
  void Ionosradar::generateInterpolationGrid(int e_nx, int e_ny, int e_nz) {
    printf("%% \tIonosradar: generateInterpolationGrid\n");

    nx = e_nx;
    ny = e_ny;
    nz = e_nz;

    unsigned int N_beams = unit_beam_dir_array.size();
    std::vector<glm::dvec3> P_box(8);

    unsigned int nd_polar = P_max_polar.size();
    double *xd_max_polar, *fd_polar;

    xd_max_polar = new double[nd_polar * 2], fd_polar = new double[nd_polar];

    printf("%% \tPoints with max range (%u) \n", nd_polar);
    for (unsigned int i = 0; i < nd_polar; i++) {

      double azm = P_max_polar[i][0], shf = 0.0;

      if (azm < 0) shf = 2.0 * M_PI;

      xd_max_polar[0 + i * 2] = azm + shf + 0.25 * M_PI;  // azmAngVal
      xd_max_polar[1 + i * 2] =
          0.5 * M_PI - P_max_polar[i][1];  // elmAngVal- -> spherical

      fd_polar[i] = P_max_polar[i][2];  // length as scalar to be interpolated

      // printf("\t\t % .5f, % .5f, % .5f \n", azm + shf, P_max_polar[i][1], P_max_polar[i][2]);
    }

    printf("%% \tComputing average beam (%lu) \n", unit_beam_dir_array.size());

    avg_dir = glm::dvec3(0.0, 0.0, 0.0);
    for (unsigned int i = 0; i < unit_beam_dir_array.size(); i++) {
      avg_dir = avg_dir + unit_beam_dir_array[i];
      if (false)
        printf("\t\t % .5f, % .5f, % .5f \n", unit_beam_dir_array[i].x,
               unit_beam_dir_array[i].y, unit_beam_dir_array[i].z);
    }

    // Normalizing
    avg_dir = (avg_dir / glm::length(avg_dir));
    printf("\tavg_dir = [% .5f,% .5f,% .5f];\n", avg_dir.x, avg_dir.y,
           avg_dir.z);

    double avg_theta = acos(avg_dir.z),
           avg_phi = atan2(avg_dir.y, avg_dir.x);  
    printf("\tavg_theta = % .5f; avg_phi = % .5f;\n\n", avg_theta, avg_phi);

    // -------------------------------------------------------------------------------
    // Computing the normalized reference as the first beam:
    printf("%% \tComputing normalized unit vectors in the K cell\n");
    double elmAngVal = 0.0, azmAngVal = 0.0,
           cx = cos(elmAngVal) * sin(azmAngVal),
           cy = cos(elmAngVal) * cos(azmAngVal), cz = sin(elmAngVal),
           dot_v0_avr;

    glm::dvec3 v0 = glm::dvec3(cx, cy, cz);
    dot_v0_avr = glm::dot(v0, avg_dir);
    glm::dvec3 vpx = v0 - dot_v0_avr * avg_dir, vpy;
    vpx = (vpx / glm::length(vpx));  // This is the reference vector for
                                     // defining the polar angle

    vpy = glm::cross(avg_dir, vpx);
    vpy = (vpy / glm::length(vpy));

    // Setting the unit reference vectors in the K cell
    K_ux = glm::dvec3(vpx.x, vpx.y, vpx.z);
    K_uy = glm::dvec3(vpy.x, vpy.y, vpy.z);
    K_uz = glm::dvec3(avg_dir.x, avg_dir.y, avg_dir.z);

    printf("\tux = [% .5f, % .5f, % .5f]; \n", K_ux.x, K_ux.y, K_ux.z);
    printf("\tuy = [% .5f, % .5f, % .5f]; \n", K_uy.x, K_uy.y, K_uy.z);
    printf("\tuz = [% .5f, % .5f, % .5f]; %% => avg_dir \n", K_uz.x, K_uz.y,
           K_uz.z);

    // Next, we project all beams into a plane that is perpendicular to the
    // avg-beam axis. This is, the resulting beams have dot(v,avg) = 1.0. Once
    // this is done, we compute the corresponding polar coords r,theta in this
    // projected plane.

    printf(
        "%% \tComputing fixed beams with range on a plane perp to avr axis "
        "(polar)\n");
    std::vector<double> cos_axis_beam_v(N_beams);
    // This is computed for building the convex hull and Fourier representation.

    vector<std::array<double, 2>> polar(
        N_beams);  // , polar_o(N_beams), transv_x(N_beams), transv_y(N_beams);

    printf("\tbeams_transv = [\t\n");
    for (unsigned int i = 0; i < N_beams; i++) {
      cos_axis_beam_v[i] = glm::dot(unit_beam_dir_array[i], avg_dir);
      double ri = 1.0 / cos_axis_beam_v[i], s_vi;

      glm::dvec3 vi = glm::dvec3(ri * unit_beam_dir_array[i].x,
                                 ri * unit_beam_dir_array[i].y,
                                 ri * unit_beam_dir_array[i].z),
                 vi_t = vi - avg_dir, vixavr;

      polar[i][0] = glm::length(vi_t);
      vixavr = glm::cross(vpx, vi_t);
      s_vi = glm::sign(glm::dot(vixavr, avg_dir));
      polar[i][1] = s_vi * acos(glm::dot(vpx, vi_t) / polar[i][0]);

      printf("\t\t % .5f, % .5f \n", polar[i][0] * cos(polar[i][1]),
             polar[i][0] * sin(polar[i][1]));
    }
    printf("\t];\n");

    printf("%% \tIonosradar::compute_quickHull_polar (x,y,ang)\n");
    vector<std::array<double, 2>> convex_hull_polar_data;

    compute_quickHull_polar(polar, convex_hull_polar_data);
    vector<std::array<double, 2>> convex_hull_data(
        convex_hull_polar_data.size());
    vector<double> convex_hull_angles(convex_hull_polar_data.size());

    double L_scale = 0.0;

    printf("\tbeams_CH = [\t\n");
    for (unsigned int i = 0; i < convex_hull_data.size(); i++) {
      double ang = convex_hull_polar_data[i][1],
             chx = convex_hull_polar_data[i][0] * cos(ang),
             chy = convex_hull_polar_data[i][0] * sin(ang);

      convex_hull_data[i][0] = chx;
      convex_hull_data[i][1] = chy;
      convex_hull_angles[i] = ang;
      printf("\t\t % .5f, % .5f\t\t% f \n", chx, chy, ang);

      L_scale += convex_hull_polar_data[i][0];
    }
    printf("\t];\n");

    L_scale = L_scale / convex_hull_data.size();

    // -------------------------------------------------------------------------------
    double F_scale = 0.9 * L_scale;
    printf("%% \tIonosradar::compute_fourier_convex_hull\n");
    compute_fourier_convex_hull(
        convex_hull_data, convex_hull_angles,
        F_scale);  // L_scale: not too tight to data points

    // Printing the coefficients
    printf("%% \tPrinting the coefficients (xA[n], xB[n], yA[n], yB[n]) \n");
    printf("\tCH_fourier_coeff = [\t\n");
    for (unsigned int n = 0; n < xA.size(); n++) {
      printf("\t\t% .10f,% .10f,% .10f, % .10f\n", xA[n], xB[n], yA[n], yB[n]);
    }
    printf("\t];\n");

    // -------------------------------------------------------------------------------
    // Optimization problem of finding the best rotation of a square to fit the
    // projected data.
    double min_angle = min_angle_sq_fit_data(convex_hull_data);
    // -------------------------------------------------------------------------------

    double Lm = 1.0e6, LM = 0.0, Rm = 1.0e6;
    RM = 0.0;

    // L is length along the average beam axis, R is perpendicular length
    // measured from L axis. We measure each scattered point and find points in
    // both faces of the outline that are perpendicular to the L axis
    // double pair_RL[2] = {0.0,0.0}; //, pair_rl[2] = {0.0,0.0}; // beams arise
    // from origin
    double cos_axis_beam, min_cos_axis_beam = 1.0;

    printf("%% \tComputing min and max radii (%lu)\n",
           unit_beam_dir_array.size());
    // Computing the bounding trapezoidal-cone: l is length and Rm, RM are
    // radius of faces.
    for (unsigned int i = 0; i < unit_beam_dir_array.size(); i++) {
      cos_axis_beam = glm::dot(unit_beam_dir_array[i], avg_dir);
      min_cos_axis_beam = min(cos_axis_beam, min_cos_axis_beam);

      double min_range = (min_max_range_array[i][0]),
             max_range = (min_max_range_array[i][1]);

      glm::dvec3 min_dir = min_range * unit_beam_dir_array[i], min_rad,
                 max_dir = max_range * unit_beam_dir_array[i], max_rad;

      double min_L = glm::dot(min_dir, avg_dir),  // min_R,
          max_L = glm::dot(max_dir, avg_dir), max_R;

      min_rad = min_dir - min_L * avg_dir;
      max_rad = max_dir - max_L * avg_dir;

      max_R = glm::length(max_rad);

      if (false)
        printf("\t\t % .5f, \t % .5f, % .5f \n", max_R,
               min_max_range_array[i][0], min_max_range_array[i][1]);

      if (min_L < Lm) Lm = min_L;
      if (max_L > LM) LM = max_L;
      if (max_R > RM) {
        RM = max_R;
      }  // pair_RL[0] = max_R; pair_RL[1] = max_L;
    }

    // similar triangles, tan_ang = RM/LM = Rm/Lm;
    Rm = Lm * RM / LM;

    // Enlarging the bounding box
    Lm = 0.90 * Lm;
    LM = 0.75 * LM;  // Less of Fourier Prism and more of Lead interpolation
    Rm = 1.00 * Rm;
    RM = 1.05 * RM;  // 1.0

    double LM_cos_beam, lm_cos_beam, R_RM,
        R_rm;  //, sq2 = sqrt(2.0); // LM_cos_beam

    LM_cos_beam = LM * min_cos_axis_beam;  // gr*LM_cos_beam+br; // gr*LM+br;
    lm_cos_beam = Lm * min_cos_axis_beam;  // gr*Lm+br;
    R_RM = RM;                             // * sq2, gr*RM+br;
    R_rm = Rm;                             // * sq2

    printf(
        "%% \tDimensions of bounding trapezoidal-cone: Lm, LM, Rm, RM, "
        "(LM-Lm)/(0.5*(Rm+RM))\n");
    printf("\tLm=% .5f; LM=% .5f; \t Rm=% .5f; RM=% .5f; C_LR=% .5f; \n", Lm,
           LM, Rm, RM, (LM - Lm) / (0.5 * (Rm + RM)));

    // double   cl = 1.5*LM/Lm; // sqrt(2.0);cv = M_PI/180.0,cr = 3.0,

    // pc   = [ R_lm, R_lm,Lm;-R_lm, R_lm,Lm;-R_lm,-R_lm,Lm; R_lm,-R_lm,Lm;
    // R_LM, R_LM,LM;-R_LM, R_LM,LM;-R_LM,-R_LM,LM; R_LM,-R_LM,LM]';

    glm::dvec3 p0 = glm::dvec3(R_rm, R_rm, lm_cos_beam),  //;
        p1 = glm::dvec3(-R_rm, R_rm, lm_cos_beam),
               p2 = glm::dvec3(-R_rm, -R_rm, lm_cos_beam),
               p3 = glm::dvec3(R_rm, -R_rm, lm_cos_beam),
               //
        p4 = glm::dvec3(R_RM, R_RM, LM_cos_beam),    // R_LM, R_LM, LM
        p5 = glm::dvec3(-R_RM, R_RM, LM_cos_beam),   //-R_LM, R_LM, LM
        p6 = glm::dvec3(-R_RM, -R_RM, LM_cos_beam),  //-R_LM,-R_LM, LM
        p7 = glm::dvec3(R_RM, -R_RM, LM_cos_beam);   // R_LM,-R_LM, LM

    // double ang_phi = 1.5*M_PI-1.0*avg_phi; //0.5*M_PI-1.0*avg_phi

    double co = cos(min_angle), so = sin(min_angle);

    P_box[0] = (p0.x * co + p0.y * so) * vpx + (-p0.x * so + p0.y * co) * vpy +
               p0.z * avg_dir;
    P_box[1] = (p1.x * co + p1.y * so) * vpx + (-p1.x * so + p1.y * co) * vpy +
               p1.z * avg_dir;
    P_box[2] = (p2.x * co + p2.y * so) * vpx + (-p2.x * so + p2.y * co) * vpy +
               p2.z * avg_dir;
    P_box[3] = (p3.x * co + p3.y * so) * vpx + (-p3.x * so + p3.y * co) * vpy +
               p3.z * avg_dir;

    P_box[4] = (p4.x * co + p4.y * so) * vpx + (-p4.x * so + p4.y * co) * vpy +
               p4.z * avg_dir;
    P_box[5] = (p5.x * co + p5.y * so) * vpx + (-p5.x * so + p5.y * co) * vpy +
               p5.z * avg_dir;
    P_box[6] = (p6.x * co + p6.y * so) * vpx + (-p6.x * so + p6.y * co) * vpy +
               p6.z * avg_dir;
    P_box[7] = (p7.x * co + p7.y * so) * vpx + (-p7.x * so + p7.y * co) * vpy +
               p7.z * avg_dir;

    printf("\tK_outline_points = [\n");
    printf("\t\t% 8.3f, % 8.3f, % 8.3f;\n", P_box[0].x, P_box[0].y, P_box[0].z);
    printf("\t\t% 8.3f, % 8.3f, % 8.3f;\n", P_box[1].x, P_box[1].y, P_box[1].z);
    printf("\t\t% 8.3f, % 8.3f, % 8.3f;\n", P_box[2].x, P_box[2].y, P_box[2].z);
    printf("\t\t% 8.3f, % 8.3f, % 8.3f;\n", P_box[3].x, P_box[3].y, P_box[3].z);
    printf("\t\t% 8.3f, % 8.3f, % 8.3f;\n", P_box[4].x, P_box[4].y, P_box[4].z);
    printf("\t\t% 8.3f, % 8.3f, % 8.3f;\n", P_box[5].x, P_box[5].y, P_box[5].z);
    printf("\t\t% 8.3f, % 8.3f, % 8.3f;\n", P_box[6].x, P_box[6].y, P_box[6].z);
    printf("\t\t% 8.3f, % 8.3f, % 8.3f;\n", P_box[7].x, P_box[7].y, P_box[7].z);
    printf("\t];\n\n");

    glm::dvec3 C_top, C_bottom, CT, CB;
    for (unsigned int i = 0; i < 4; i++) {
      C_bottom = C_bottom + P_box[i];
      C_top = C_top + P_box[i + 4];
    }

    CB = P_box[0] - 0.25 * C_bottom;
    CT = P_box[4] - 0.25 * C_top;

    double Rb = glm::length(CB);  //, Rt = glm::length(CT);

    // double x,y,z;
    // unsigned int i;

    string spoints = "scattered_points";  // output_folder
    output_radar_data_vtk(spoints, P_data);

    double LC = LM_cos_beam,                 // LM*(1.0-cos(max_angle)),
        RC = 1.0 * sqrt(LM * LM - LC * LC);  // LM*sin(max_angle);
    std::vector<glm::dvec3> Pc_box(8);
    for (unsigned int i = 0; i < 4; i++) {
      double pz = glm::dot(P_box[i],
                           avg_dir);  // , aux = pz/( glm::length(P_box[i])*1.0
                                      // ), angle = acos ( aux );
      glm::dvec3 Pp = pz * avg_dir, Pt = P_box[i] - Pp;

      Pt = Pt * (1.0 / glm::length(Pt));

      // Rm = Lm*RM/LM;
      Pc_box[i] = Rb * Pt + Lm * avg_dir;
      Pc_box[i + 4] = RC * Pt + LC * avg_dir;
    }

    for (unsigned int i = 0; i < 4; i++) {
      double px = glm::dot(P_box[i], K_ux), py = glm::dot(P_box[i], K_uy),
             pz = glm::dot(P_box[i], K_uz), oc = atan2(py, px),
             fr = fourier_radius(oc);

      glm::dvec3 Pp = pz * avg_dir, Pt = P_box[i] - Pp;

      Pt = Pt * (1.0 / glm::length(Pt));

      // Rm = Lm*RM/LM;
      Pc_box[i] = Rb * fr * Pt + Lm * avg_dir;
      Pc_box[i + 4] = RC * fr * Pt + LC * avg_dir;
    }

    // We fix the top layer to match the Shepard interpolant of max range data
    double* pbR;
    pbR = new double[1];
    pbR[0] = LM;

    unsigned int n_poly = 3 + 1;
    double dl = 1.0 / (n_poly - 1.0);

    std::vector<double> top_poly_points(n_poly),
        top_poly_coeff(n_poly * n_poly);

    for (unsigned int i = 0; i < n_poly; i++) {
      double xip = i * dl;
      top_poly_points[i] = xip;
      for (unsigned int j = 0; j < n_poly; j++) {
        double eta = j * dl;

        top_poly_coeff[i * n_poly + j] = top_lead_interpolation(
            xip, eta, Pc_box, nd_polar, xd_max_polar, fd_polar);
      }
    }

    scaled3D_TI_Fourier_CH_Poly_top_grid(nx, ny, nz, Pc_box, n_poly,
                                         top_poly_points, top_poly_coeff, pi,
                                         Rb, RC, LM);

    ni = pi.size();
    xi = new double[3 * ni];

    for (unsigned int i = 0; i < ni; i++) {
      xi[0 + 3 * i] = pi[i][0];
      xi[1 + 3 * i] = pi[i][1];
      xi[2 + 3 * i] = pi[i][2];
    }

    P_box.clear();
    Pc_box.clear();

    delete[] xd_max_polar;
    delete[] fd_polar;
  }
  // -----------------------------------------------------------------------------------
// Shepard interpolation for Ion data
void Ionosradar::sparseInterpolationAlgorithm_large (
      int timestep, std::vector<std::array<double, 8>>& data_array,
      string fname) {
    R = 1000.0;
    //double c_nu = smoothness;
    double p = shepard_parameter;
    string str_der;
    const unsigned NS = 5U;
    
    // neli, nOpi, vsi, Tii, Tei
    std::vector<std::array<double, NS>> data_interp(ni);
	
    {  // Vectorized
      unsigned ND = data_array.size();
      constexpr unsigned NS = 5U;
      unsigned DIMS = 3U;
      unsigned VEC_WIDTH = 8U;
      unsigned NX = nx;
      unsigned NY = ny;
      unsigned NZ = nz;
      unsigned NI = ni;
      float* _data = new float[ND * (DIMS + NS)];
      float* _xi = new float[3U * NI];
      float* Fi = allocateAlignedResult(NI, NS, VEC_WIDTH);
      unsigned NC = NUM_THREADS; //8; //omp_get_num_threads();

      for (unsigned i = 0; i < NI; i++) {
        _xi[0 * NI + i] = pi[i][0];
        _xi[1 * NI + i] = pi[i][1];
        _xi[2 * NI + i] = pi[i][2];
      }

      for (unsigned i = 0; i < ND; i++) {
        for (unsigned j = 0; j < (DIMS + NS); j++) {
          _data[i * (DIMS + NS) + j] = data_array[i][j];
        }
      }

      spInterpMultiVecParallel(ND, NS, NX, NY, NZ, NI, _xi, _data, NC, Fi);
      unpackInterpolation(NI, NS, VEC_WIDTH, Fi, data_interp);
    }
    // ***************************************************************************

    // Printing routine: Verbose with -v option
    for (unsigned int i = 0; i < ni; i++) {  // ni

      double neli, nOpi, vsi, Tii, Tei;

      neli = data_interp[i][0];
      nOpi = data_interp[i][1];
      vsi  = data_interp[i][2];
      Tii  = data_interp[i][3];
      Tei  = data_interp[i][4];

      if (verbose)
        printf("i=%d:\tISR(%5.3f,\t%5.3f,\t%.3e,\t%.3e,\t%.3e)\n", i, neli,
               nOpi, vsi, Tii, Tei);
    }

    string sname = "_p" + std::to_string((int)p);
    write_interpolated_time_selected_vals_vtk("SEL" + sname, num_scalar_data,
                                              timestep, data_interp);

    data_interp.clear();
    // printf("Time %d NEL measured: %.3f seconds.\n", timestep, t_interp);
  }  
  // -----------------------------------------------------------------------------------
//********************** single scalar version **********************
  double* Ionosradar::shepard_interp_der_2d(int p, unsigned int nd, double* xd, double* fd,
                                unsigned int ni, double* xi,
                                bool estimate_derivatives, double R, double RM,
                                double c_nu)

  //*******************************************************************
  //
  //  Purpose:
  //    SHEPARD_INTERP_ND evaluates a multidimensional Shepard interpolant.
  //
  //  Licensing:
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //    22 June 2022
  //
  //  Author:
  //    John Burkardt, Juan Carlos Araújo
  //
  //  Reference:
  //
  //    Donald Shepard,
  //    A two-dimensional interpolation function for irregularly spaced data,
  //    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
  //    ACM, pages 517-524, 1969.
  //
  //  Parameters:
  //
  //    Input,  int M, the spatial dimension.
  //    Input,  int ND, the number of data points.
  //    Input,  double XD[M*ND], the data points.
  //    Input,  double fd[ND], the data values.
  //    Input,  int NI, the number of interpolation points.
  //    Input,  double XI[M*NI], the interpolation points.
  //    Output, double SHEPARD_INTERP_ND[NI], the interpolated values.
  //
  //	  R, interpolate is point is to a distance R
  //    RM, threshold, if |xi| > RM, then we ...
  //
  //************************************************************************************
  {
    
    double s;
    double t;
    double *w, *Fx, *Fy, *di;  // *Fz,
    int z;
    double* Fi;

    unsigned int dim = 2;

    Fi = new double[ni];
    w = new double[nd];
    di = new double[nd];

    Fx = new double[nd];
    Fy = new double[nd];

    double nu = 0.0, aux;

    // Correction by including derivative estimation
    double max_Fx2 = 0.0, max_Fy2 = 0.0, min_f = 1.0e20, max_f = -1.0e20;

    if (estimate_derivatives) {
      double WxdF_t2, t2, wij, max_DF;

      for (unsigned int ii = 0; ii < nd; ii++) {
        if (!isnan(fd[ii])) {
          s = 0.0;

          Fx[ii] = 0.0;
          Fy[ii] = 0.0;

          for (unsigned int j = 0; j < nd; j++) {
            if (!isnan(fd[j])) {
              min_f = min(min_f, fd[j]);
              max_f = max(max_f, fd[j]);
              if (j != ii) {
                t2 = 0.0;
                for (unsigned int i2 = 0; i2 < dim; i2++) {
                  t2 = t2 + pow(xd[i2 + j * dim] - xd[i2 + ii * dim],
                                2);  // xd[i2+j*dim] - xd[i2+ii*dim]
                }
                wij = 1.0 / t2;
                t = sqrt(t2);
                s = s + wij;

                WxdF_t2 = wij * (fd[j] - fd[ii]) / t2;

                Fx[ii] =
                    Fx[ii] +
                    WxdF_t2 * (xd[0 + j * dim] -
                               xd[0 + ii * dim]);  // xd[0+j*dim] - xd[0+ii*dim]
                Fy[ii] =
                    Fy[ii] +
                    WxdF_t2 * (xd[1 + j * dim] -
                               xd[1 + ii * dim]);  // xd[1+j*dim] - xd[1+ii*dim]
              }  // j != ii
            }  // if not nan
          }  // j
          if (s > 1e-14) {
            // printf ("** s **: i=%u, ii=%u\n",i,j);

            Fx[ii] = Fx[ii] / s;
            Fy[ii] = Fy[ii] / s;
          } else {
            Fx[ii] = 0.0;
            Fy[ii] = 0.0;
          }
          max_Fx2 = max(max_Fx2, Fx[ii] * Fx[ii]);
          max_Fy2 = max(max_Fy2, Fy[ii] * Fy[ii]);  // fixed
        }  // if not nan
      }  // ii

      if (max_Fx2 + max_Fy2 > 0) {
        max_DF = sqrt(max_Fx2 + max_Fy2);
        nu = c_nu * (max_f - min_f) / max_DF;
      } else {
        nu = 0.0;
      }
    }  // if derivatives

    for (unsigned int i = 0; i < ni; i++) {
      t = 0.0;
      // bool inside_RM = true;

      z = -1;
      for (unsigned int j = 0; j < nd; j++) {
        min_f = min(min_f, fd[j]);
        max_f = max(max_f, fd[j]);

        t = 0.0;

        for (unsigned int i2 = 0; i2 < dim; i2++) {
          t = t + pow(xi[i2 + i * dim] - xd[i2 + j * dim],
                      2);  // xi[i2+i*dim] - xd[i2+j*dim]
        }
        di[j] = sqrt(t);
        w[j] = di[j];

        if (w[j] == 0.0) {
          z = j;
          break;
        }
      }  // j

      if (z != -1) {  // singular case where di[j] = 0.0
        for (unsigned int j = 0; j < nd; j++) {
          w[j] = 0.0;
        }
        w[z] = 1.0;
      } else {  // case where di[j] > 0.0
        s = 0.0;
        for (unsigned int j = 0; j < nd; j++) {
          if (isnan(fd[j])) {  // missing data from ISR
            w[j] = 0.0;
          } else {
            
            aux = 1. / w[j];

            w[j] = pow(aux, p);  
            s = s + w[j];
          }
        }  // j

        for (unsigned int j = 0; j < nd; j++) {
          w[j] = w[j] / s;
        }
      }

      Fi[i] = 0.0;  // 0.5*(min_f+max_f);
      for (unsigned int j = 0; j < nd; j++) {
        if (!isnan(fd[j])) Fi[i] += w[j] * fd[j];
      }
      if (estimate_derivatives) {
        // Fi_new(i) = Fi(i);
        for (unsigned int j = 0; j < nd; j++) {  // fd[j][i2] - fd[ii][i2]
          Fi[i] =
              Fi[i] +
              w[j] *
                  (Fx[j] * (xi[0 + i * dim] -
                            xd[0 + j * dim]) +  //(xi(i,1) - xd(j,1)) + //
                                                // xd[0+j*dim]
                   Fy[j] *
                       (xi[1 + i * dim] -
                        xd[1 + j * dim])  //(xi(i,2) - xd(j,2)) + // xd[1+j*dim]
                   ) *
                  nu / (nu + di[j]);
        }
      }  // if derivatives

      if (isnan(Fi[i])) {
        printf("\t\t% .3f,% .3f,% .5f \t\t %% shepard_interp_der_2d\n",
               xi[0 + i * dim], xi[1 + i * dim], Fi[i]);
      }

      if (abs(Fi[i]) <
          1e-14) {  // in some cases the radar data is missing in the
                    // first data of the beam ... then we wrongly
                    // assign Fi = 0. We should assign min_f
        Fi[i] = min_f;
      }
    }  // xi

    delete[] w;
    delete[] di;
    delete[] Fx;
    delete[] Fy;

    return Fi;
  }
 // -------------------------------------------------------------------------------------
double Ionosradar::fourier_radius(double t) {
    unsigned int Nf = xA.size() - 1;  //, count = 1; // xA.size()-1;
    double xf = 0.0, yf = 0.0, rf, snt, cnt;

    xf = xA[0];
    yf = yA[0];
    for (unsigned int n = 1; n <= Nf; n++) {  // remove the last entry
      snt = sin(n * t);
      cnt = cos(n * t);
      xf = xf + xA[n] * cnt + xB[n] * snt;
      yf = yf + yA[n] * cnt + yB[n] * snt;
    }

    rf = sqrt(xf * xf + yf * yf);

    return rf;
  }

double* Ionosradar::fourier_radius_xy(double t) {
    unsigned int Nf = xA.size() - 1;      // xA.size()-1;
    double xf = 0.0, yf = 0.0, snt, cnt;  // r

    xf = xA[0];
    yf = yA[0];
    for (unsigned int n = 1; n <= Nf; n++) {  // remove the last entry
      snt = sin(n * t);
      cnt = cos(n * t);
      xf = xf + xA[n] * cnt + xB[n] * snt;
      yf = yf + yA[n] * cnt + yB[n] * snt;
    }
    
    static double v[2];

    v[0] = xf;
    v[1] = yf;

    return v;
  }

glm::dvec3 Ionosradar::fourier_xy(double t) {
    unsigned int Nf = xA.size() - 1;      // xA.size()-1;
    double xf = 0.0, yf = 0.0, snt, cnt;  // r

    xf = xA[0];
    yf = yA[0];
    for (unsigned int n = 1; n <= Nf; n++) {  // remove the last entry
      snt = sin(n * t);
      cnt = cos(n * t);
      xf = xf + xA[n] * cnt + xB[n] * snt;
      yf = yf + yA[n] * cnt + yB[n] * snt;
    }
    // r = sqrt (xf*xf + yf*yf);
    glm::dvec3 v(xf, yf, 0.0);

    return v;
  }
    
// -----------------------------------------------------------------------------------
  // std::array<double,3>
  double Ionosradar::top_lead_interpolation(double x, double y, std::vector<glm::dvec3>& V,
                                unsigned int nd_e, double* xd_e,
                                double* fd_e) {  // Shepard interpolation data

    double phi_s[4], ang_xy, ang_z_plane, r_xy_plane;
    unsigned int idx[4];

    glm::dvec3 Xf;  // X, X0, P

    // mapping only to the upper face
    phi_s[0] = (1.0 - x) * (1.0 - y);  // v4
    phi_s[1] = x * (1.0 - y);          // v5
    phi_s[2] = x * y;                  // v6
    phi_s[3] = (1.0 - x) * y;          // v7

    idx[0] = 4;
    idx[1] = 5;
    idx[2] = 6;
    idx[3] = 7;  //

    Xf.x = 0.0;
    Xf.y = 0.0;
    Xf.z = 0.0;
    for (unsigned int i = 0; i < 4; i++) {
      Xf = Xf + phi_s[i] * V[idx[i]];
    }

    glm::dvec3 Pf = Xf;  // P,,,,

    double* R;
    R = new double[1];
    R[0] = 1.0;  // = LM;

    r_xy_plane = sqrt(Pf.x * Pf.x + Pf.y * Pf.y);
    ang_xy = atan2(Pf.y, Pf.x);
    ang_z_plane = atan2(r_xy_plane, Pf.z);

    double azm = ang_xy, shf = 0.0, Rf;
    if (azm < 0) shf = 2.0 * M_PI;

    double* xit;
    xit = new double[2 * 1];
    xit[0] = ang_xy + shf;
    xit[1] = ang_z_plane;
    
    R = shepard_interp_der_2d(2.0, nd_e, xd_e, fd_e, 1, xit, true, 1000.0,
                              1000.0,
                              100 * smoothness);  // shepard_parameter
    Rf = R[0];
    
    return Rf;  // point;
  }

// -----------------------------------------------------------------------------------
  std::array<double, 3> Ionosradar::transfinite_interpolation_Fourier_CH_Poly_top(
      double x, double y, double z, std::vector<glm::dvec3>& V,
      unsigned int n_poly, std::vector<double>& top_poly_points,
      std::vector<double>& top_poly_coeff, double R_rm, double R_RM,
      double LM) {
    // here, x,y,z belong to the unit cube [0,1]^3 and we map to the physical
    // brick element
    double phi[8], phi_s[4], phi_sx0[4], phi_sx1[4], phi_sy0[4], phi_sy1[4], xc,
        yc, zc, oc, rho, ang_xy, ang_z_plane,
        r_xy_plane;  // rf, rc = 1.0, rho_n, theta
    unsigned int idx[4];

    glm::dvec3 uv, rho_v, aux_v;

    phi[0] = (1.0 - x) * (1.0 - y) * (1.0 - z);
    phi[1] = x * (1.0 - y) * (1.0 - z);
    phi[2] = x * y * (1.0 - z);
    phi[3] = (1.0 - x) * y * (1.0 - z);

    // mapping only to the upper face
    phi_s[0] = (1.0 - x) * (1.0 - y);  // v4
    phi_s[1] = x * (1.0 - y);          // v5
    phi_s[2] = x * y;                  // v6
    phi_s[3] = (1.0 - x) * y;          // v7

    phi[4] = phi_s[0] * z;
    phi[5] = phi_s[1] * z;
    phi[6] = phi_s[2] * z;
    phi[7] = phi_s[3] * z;

    glm::dvec3 X, Xf, X0, P;

    // Physical brick
    X.x = 0.0;
    X.y = 0.0;
    X.z = 0.0;
    for (unsigned int i = 0; i < 8; i++) {
      X = X + phi[i] * V[i];
    }

    P = 1.0 * X;

    // Top face z=1
    // ------------------------------------------------------------------
    idx[0] = 4;
    idx[1] = 5;
    idx[2] = 6;
    idx[3] = 7;  // {};

    Xf.x = 0.0;
    Xf.y = 0.0;
    Xf.z = 0.0;
    for (unsigned int i = 0; i < 4; i++) {
      Xf = Xf + phi_s[i] * V[idx[i]];
    }

    glm::dvec3 Pf = Xf;  // P,,,,

    // double* R; R = new double[1]; R[0] = 1.0; // = LM;

    r_xy_plane = sqrt(Pf.x * Pf.x + Pf.y * Pf.y);
    ang_xy = atan2(Pf.y, Pf.x);
    ang_z_plane = atan2(r_xy_plane, Pf.z);

    // double azm = ang_xy; //, shf = 0.0;
    // if (azm < 0) shf = 2.0*M_PI;

    double Rf = 0.0, Rp;
    for (unsigned int i = 0; i < n_poly; i++) {
      for (unsigned int j = 0; j < n_poly; j++) {
        unsigned int k = i * n_poly + j;
        Rf = Rf + top_poly_coeff[k] * Lshape(i, top_poly_points, x) *
                      Lshape(j, top_poly_points, y);  // p;
      }
    }

    Rp = 1.05 * Rf;

    X0.x = Rp * sin(ang_z_plane) * cos(ang_xy) - Xf.x;
    X0.y = Rp * sin(ang_z_plane) * sin(ang_xy) - Xf.y;
    X0.z = Rp * cos(ang_z_plane) - Xf.z;

    P = P + X0 * z;

    // // mapping only to the face x=0
    // -----------------------------------------------
    idx[0] = 0;
    idx[1] = 3;
    idx[2] = 4;
    idx[3] = 7;  // {0,3,4,7}; -> 1,4,5,8

    phi_sx0[0] = (1.0 - y) * (1.0 - z);
    phi_sx0[1] = y * (1.0 - z);
    phi_sx0[2] = (1.0 - y) * z;
    phi_sx0[3] = y * z;

    Xf.x = 0.0;
    Xf.y = 0.0;
    Xf.z = 0.0;
    for (unsigned int i = 0; i < 4; i++) {
      Xf = Xf + phi_sx0[i] * V[idx[i]];
    }

    xc = glm::dot(Xf, K_ux);
    yc = glm::dot(Xf, K_uy);
    zc = glm::dot(Xf, K_uz);

    uv = zc * K_uz;  // vec along the axis
    rho = R_rm +
          (R_RM - R_rm) * z;  // conical: linear interpolation between the radii

    oc = atan2(yc, xc);
    aux_v = rho * fourier_xy(oc);
    rho_v = aux_v.x * K_ux + aux_v.y * K_uy;

    X0.x = rho_v.x + uv.x - Xf.x;
    X0.y = rho_v.y + uv.y - Xf.y;
    X0.z = rho_v.z + uv.z - Xf.z;

    P = P + X0 * (1.0 - x);

    // // mapping only to the face x=1
    // -----------------------------------------------
    idx[0] = 1;
    idx[1] = 2;
    idx[2] = 5;
    idx[3] = 6;  // {1,2,5,6}; -> 2,3,6,7

    phi_sx1[0] = (1.0 - y) * (1.0 - z);
    phi_sx1[1] = y * (1.0 - z);
    phi_sx1[2] = (1.0 - y) * z;
    phi_sx1[3] = y * z;

    Xf.x = 0.0;
    Xf.y = 0.0;
    Xf.z = 0.0;
    for (unsigned int i = 0; i < 4; i++) {
      Xf = Xf + phi_sx1[i] * V[idx[i]];
    }

    xc = glm::dot(Xf, K_ux);
    yc = glm::dot(Xf, K_uy);
    zc = glm::dot(Xf, K_uz);

    uv = zc * K_uz;  // vec along the axis
    rho = R_rm +
          (R_RM - R_rm) * z;  // conical: linear interpolation between the radii

    oc = atan2(yc, xc);
    aux_v = rho * fourier_xy(oc);
    rho_v = aux_v.x * K_ux + aux_v.y * K_uy;

    X0.x = rho_v.x + uv.x - Xf.x;
    X0.y = rho_v.y + uv.y - Xf.y;
    X0.z = rho_v.z + uv.z - Xf.z;

    P = P + X0 * x;

    // // mapping only to the face y=0
    // -----------------------------------------------
    idx[0] = 0;
    idx[1] = 1;
    idx[2] = 4;
    idx[3] = 5;  // {0,1,4,5}; -> 1,2,5,6
    phi_sy0[0] = (1.0 - x) * (1.0 - z);
    phi_sy0[1] = x * (1.0 - z);
    phi_sy0[2] = (1.0 - x) * z;
    phi_sy0[3] = x * z;

    Xf.x = 0.0;
    Xf.y = 0.0;
    Xf.z = 0.0;
    for (unsigned int i = 0; i < 4; i++) {
      Xf = Xf + phi_sy0[i] * V[idx[i]];
    }

    xc = glm::dot(Xf, K_ux);
    yc = glm::dot(Xf, K_uy);
    zc = glm::dot(Xf, K_uz);

    uv = zc * K_uz;  // vec along the axis
    // rho_v = Xf - uv;
    rho = R_rm +
          (R_RM - R_rm) * z;  // conical: linear interpolation between the radii

    oc = atan2(yc, xc);
    aux_v = rho * fourier_xy(oc);
    rho_v = aux_v.x * K_ux + aux_v.y * K_uy;

    X0.x = rho_v.x + uv.x - Xf.x;
    X0.y = rho_v.y + uv.y - Xf.y;
    X0.z = rho_v.z + uv.z - Xf.z;

    P = P + X0 * (1.0 - y);

    // // mapping only to the face y=1
    // -----------------------------------------------
    idx[0] = 2;
    idx[1] = 3;
    idx[2] = 6;
    idx[3] = 7;  // {2,3,6,7}; -> 3,4,7,8
    phi_sy1[0] = x * (1.0 - z);
    phi_sy1[1] = (1.0 - x) * (1.0 - z);
    phi_sy1[2] = x * z;
    phi_sy1[3] = (1.0 - x) * z;

    Xf.x = 0.0;
    Xf.y = 0.0;
    Xf.z = 0.0;
    for (unsigned int i = 0; i < 4; i++) {
      Xf = Xf + phi_sy1[i] * V[idx[i]];
    }

    xc = glm::dot(Xf, K_ux);
    yc = glm::dot(Xf, K_uy);
    zc = glm::dot(Xf, K_uz);

    uv = zc * K_uz;  // vec along the axis
    // rho_v = Xf - uv;
    rho = R_rm +
          (R_RM - R_rm) * z;  // conical: linear interpolation between the radii

    oc = atan2(yc, xc);
    aux_v = rho * fourier_xy(oc);
    rho_v = aux_v.x * K_ux + aux_v.y * K_uy;

    X0.x = rho_v.x + uv.x - Xf.x;
    X0.y = rho_v.y + uv.y - Xf.y;
    X0.z = rho_v.z + uv.z - Xf.z;

    P = P + X0 * y;

    std::array<double, 3> point = {P.x, P.y, P.z};
    return point;
  }
  // -------------------------------------------------------------------------------------
void Ionosradar::scaled3D_TI_Fourier_CH_Poly_top_grid (
      unsigned int nx, unsigned int ny, unsigned int nz,
      std::vector<glm::dvec3>& V,
      // unsigned int nd_e, double *xd_e, double *fd_e,
      unsigned int n_poly, std::vector<double>& top_poly_points,
      std::vector<double>& top_poly_coeff,
      std::vector<std::array<double, 3>>& P, double R_rm, double R_RM,
      double LM) {
    double dx = 1.0 / (nx - 1);
    double dy = 1.0 / (ny - 1);
    double dz = 1.0 / (nz - 1);
    //,xm,ym,zm,sx,sy,sz, pp = 1.25;
    std::array<double, 3> point;

    for (unsigned int i = 0; i < nx; i++) {
      for (unsigned int j = 0; j < ny; j++) {
        for (unsigned int k = 0; k < nz; k++) {
          point = transfinite_interpolation_Fourier_CH_Poly_top(
              dx * i, dy * j, dz * k, V, n_poly, top_poly_points,
              top_poly_coeff, R_rm, R_RM, LM);

          P.push_back(point);
        }
      }
    }
  }
  // -------------------------------------------------------------------------------------
  void Ionosradar::collect_BeamData (
      H5::Group& beamGroup, std::vector<std::array<double, 3>>& points_array,
      std::vector<std::array<double, 3>>& points_max_array,
      std::vector<std::array<double, 3>>& points_max_array_polar,
      std::vector<glm::dvec3>& unit_beam_dir_array,
      std::vector<std::array<double, 2>>& min_max_range_array,
      std::array<double, 6>& beam_bounding_box) {
    H5::Group params1d = beamGroup.openGroup(PAR_1D);
    H5::Group params2d = beamGroup.openGroup(PAR_2D);
    H5::DataSet range = beamGroup.openDataSet(RANGE_DATA);
    H5::DataSet time = beamGroup.openDataSet(TIME_DATA);

    H5::DataSet nel_beamDataSet =
        params2d.openDataSet(BEAM_DATA_NEL);  // BEAM_DATA_UNEL

    H5::DataSet azmAng = params1d.openDataSet(ANG_AZM_DATA);
    H5::DataSet elmAng = params1d.openDataSet(ANG_ELM_DATA);

    hsize_t numOfRange = getDimensions(range).cols;
    hsize_t numOfTime = getDimensions(time).cols;

    hsize_t numOfazmAng = getDimensions(azmAng).cols;
    hsize_t numOfelmAng = getDimensions(elmAng).cols;

    Dim numOfDimSamples = getDimensions(nel_beamDataSet);

    double* azmAngData = read1DData<double>(H5::PredType::NATIVE_DOUBLE, azmAng, 0,
                                            numOfazmAng, numOfazmAng);
    double* elmAngData = read1DData<double>(H5::PredType::NATIVE_DOUBLE, elmAng, 0,
                                            numOfelmAng, numOfelmAng);
    double* rangeData = read1DData<double>(H5::PredType::NATIVE_DOUBLE, range, 0,
                                           numOfRange, numOfRange);

    int* timeData =
        read1DData<int>(H5::PredType::NATIVE_INT, time, 0, numOfTime, numOfTime);

    // double lastLen = rangeData[numOfRange - 1];

    time_measurements = numOfDimSamples.rows;

    H5::DataSet beamid_set = params1d.openDataSet("beamid");
    hsize_t numOfBeamsData = getDimensions(beamid_set).cols;
    int* beamid = read1DData<int>(H5::PredType::NATIVE_INT, beamid_set, 0,
                                  numOfBeamsData, numOfBeamsData);

    double radians = M_PI / 180.0;

    double m_x = 1e8, M_x = -1e8, m_y = 1e8, M_y = -1e8, m_z = 1e8, M_z = -1e8;

    min_alt = m_z;
    max_alt = M_z;

    // printf (" *********** angles, numOfRange = %d ************** \n",
    // (int)numOfRange);

    {
      int j = 0;
      // printf("%%**** j=%d\n", j);
      //  double radians = M_PI/180.0;
      double elmAngVal = radians * (elmAngData[j]);
      double azmAngVal = radians * (azmAngData[j]);

      // printf (" \t\t\t elmAngVal = %.3f, azmAngVal = %.3f \n ",
      // elmAngData[j], azmAngData[j] );

      double cx = cos(elmAngVal) * sin(azmAngVal);
      double cy = cos(elmAngVal) * cos(azmAngVal);
      double cz = sin(elmAngVal);

      // glm::dvec3 start = { lastLen,elmAngVal,azmAngVal };
      glm::dvec3 end = glm::dvec3(cx, cy, cz), z_axis = glm::dvec3(0, 0, 1);
      unit_beam_dir_array.push_back(end);

      // closest axis to z_axis
      double local_res_z = abs(glm::dot(end, z_axis));

      if (local_res_z > res_z) {
        res_z = local_res_z;
        beam_closest_z_axis = beamGroup;
        unit_beam_test_array = end;
        // printf (" +++++++++ res_z=%.5f +++++++++++", res_z);
      }

      // unsigned int id, index0, indexf, num_range;
      // double elm, azm;

      // Storing beam info
      BeamInfo current_beam;

      current_beam.id = beamid[0];
      current_beam.num_range = numOfRange;
      current_beam.elm = elmAngData[j];
      current_beam.azm = azmAngData[j];

      // current_beam. = ;
      // beam_info belongs to Ionosradar
      beam_info.push_back(current_beam);

      std::array<double, 2> min_max_beam = {rangeData[0],
                                            rangeData[numOfRange - 1]};
      min_max_range_array.push_back(min_max_beam);

      for (unsigned int i = 0; i < numOfRange; i++) {  // range
        double length = rangeData[i], min_range = 1.0e6, max_range = 0.0;

        if (length < min_range) min_range = length;
        if (length > max_range) max_range = length;

        glm::dvec3 P = length * end;

        m_x = min(m_x, P.x);
        M_x = max(M_x, P.x);
        m_y = min(m_y, P.y);
        M_y = max(M_y, P.y);
        m_z = min(m_z, P.z);
        M_z = max(M_z, P.z);

        // P_data.push_back( {P.x,P.y,P.z} );
        points_array.push_back({P.x, P.y, P.z});

        if (i == numOfRange - 1) {
          points_max_array.push_back({P.x, P.y, P.z});
          points_max_array_polar.push_back({azmAngVal, elmAngVal, length});
        }
        min_max_beam[0] = min_range;
        min_max_beam[1] = max_range;
      }

      beam_bounding_box[0] = m_x;
      beam_bounding_box[1] = M_x;
      beam_bounding_box[2] = m_y;
      beam_bounding_box[3] = M_y;
      beam_bounding_box[4] = m_z;
      beam_bounding_box[5] = M_z;

      min_alt = min(min_alt, m_z);
      max_alt = max(max_alt, M_z);
    }

    delete[] azmAngData;
    delete[] elmAngData;
    delete[] rangeData;
    delete[] timeData;

  }  // collect_BeamData
  // -------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------
  void Ionosradar::collect_BeamData_Large(H5::Group& beamGroup,
                              std::vector<std::array<double, 8>>& data_array,
                              unsigned int time_step) {
    H5::Group params1d = beamGroup.openGroup(PAR_1D);
    H5::Group params2d = beamGroup.openGroup(PAR_2D);
    H5::DataSet range = beamGroup.openDataSet(RANGE_DATA);
    H5::DataSet time = beamGroup.openDataSet(TIME_DATA);

    H5::DataSet azmAng = params1d.openDataSet(ANG_AZM_DATA);
    H5::DataSet elmAng = params1d.openDataSet(ANG_ELM_DATA);

    hsize_t numOfRange = getDimensions(range).cols;
    hsize_t numOfTime = getDimensions(time).cols;

    hsize_t numOfazmAng = getDimensions(azmAng).cols;
    hsize_t numOfelmAng = getDimensions(elmAng).cols;

    H5::DataSet nel_beamDataSet =
        params2d.openDataSet(BEAM_DATA_NEL);  // BEAM_DATA_UNEL
    Dim numOfDimSamples = getDimensions(nel_beamDataSet);

    double* azmAngData = read1DData<double>(H5::PredType::NATIVE_DOUBLE, azmAng, 0,
                                            numOfazmAng, numOfazmAng);
    double* elmAngData = read1DData<double>(H5::PredType::NATIVE_DOUBLE, elmAng, 0,
                                            numOfelmAng, numOfelmAng);
    double* rangeData = read1DData<double>(H5::PredType::NATIVE_DOUBLE, range, 0,
                                           numOfRange, numOfRange);

    int* timeData =
        read1DData<int>(H5::PredType::NATIVE_INT, time, 0, numOfTime, numOfTime);

    H5::DataSet pop_beamDataSet, vo_beamDataSet, Ti_beamDataSet, Te_beamDataSet;
    Array2D<double> nel_data, pop_data, vo_data, ti_data, te_data;

    nel_data = read2DData<double>(H5::PredType::NATIVE_DOUBLE, nel_beamDataSet,
                                  {0, 0}, numOfDimSamples, numOfDimSamples);

    if (num_scalar_data > 4) {
      pop_beamDataSet = params2d.openDataSet(BEAM_DATA_POp);
      vo_beamDataSet = params2d.openDataSet(BEAM_DATA_VO);
      Ti_beamDataSet = params2d.openDataSet(BEAM_DATA_TI);
      Te_beamDataSet = params2d.openDataSet(BEAM_DATA_TE);

      pop_data = read2DData<double>(H5::PredType::NATIVE_DOUBLE, pop_beamDataSet,
                                    {0, 0}, numOfDimSamples, numOfDimSamples);
      vo_data = read2DData<double>(H5::PredType::NATIVE_DOUBLE, vo_beamDataSet,
                                   {0, 0}, numOfDimSamples, numOfDimSamples);
      ti_data = read2DData<double>(H5::PredType::NATIVE_DOUBLE, Ti_beamDataSet,
                                   {0, 0}, numOfDimSamples, numOfDimSamples);
      te_data = read2DData<double>(H5::PredType::NATIVE_DOUBLE, Te_beamDataSet,
                                   {0, 0}, numOfDimSamples, numOfDimSamples);
    }

    // double lastLen = rangeData[numOfRange - 1];

    double radians = M_PI / 180.0;
    {
      int j = time_step;

      double elmAngVal = radians * (elmAngData[j]);
      double azmAngVal = radians * (azmAngData[j]);

      double cx = cos(elmAngVal) * sin(azmAngVal);
      double cy = cos(elmAngVal) * cos(azmAngVal);
      double cz = sin(elmAngVal);

      glm::dvec3 end = glm::dvec3(cx, cy, cz);
      unit_beam_dir_array.push_back(end);

      for (unsigned int i = 0; i < numOfRange; i++) {
        double length = rangeData[i];

        glm::dvec3 P = length * end;

        // Here, we should allow nan numbers and in the interpolation scheme, we
        // apply the procedure for missing data.

        {
          double ne = std::numeric_limits<double>::quiet_NaN(),
                 nOp = std::numeric_limits<double>::quiet_NaN(),
                 Vs = std::numeric_limits<double>::quiet_NaN(),
                 Ti = std::numeric_limits<double>::quiet_NaN(),
                 Te = std::numeric_limits<double>::quiet_NaN();

          ne = nel_data[i][j];
          if (num_scalar_data > 4) {
            nOp = pop_data[i][j];
            Vs = vo_data[i][j];
            Ti = ti_data[i][j];
            Te = te_data[i][j];
          }
          data_array.push_back({P.x, P.y, P.z, ne, nOp, Vs, Ti, Te});
        }
      }
    }

    delete[] azmAngData;
    delete[] elmAngData;
    delete[] rangeData;
    delete[] timeData;

    nel_data.clear();
    if (num_scalar_data > 4) {
      pop_data.clear();
      vo_data.clear();
      ti_data.clear();
      te_data.clear();
    }
  }
  // -----------------------------------------------------------------------------------
   
float Ionosradar::ratioNANne(const std::vector<std::array<double, 8>>& radar_data) {
  float ratio = 0.0f;
  for (unsigned i = 0U; i < radar_data.size(); i++) {
    if (std::isnan(radar_data[i][3])) ratio += 1.0f;
  }
  return ratio / static_cast<float>(radar_data.size());
}

void Ionosradar::produceNANResult(int timestep) {
  constexpr unsigned NS = 5U;
  std::vector<std::array<double, NS>> data_interp(ni);

  // Printing routine: Verbose with -v option. Plus, setting the array to nan.
  for (unsigned i = 0U; i < ni; i++) {
    data_interp[i] = {std::numeric_limits<double>::quiet_NaN(), 
                      std::numeric_limits<double>::quiet_NaN(), 
                      std::numeric_limits<double>::quiet_NaN(), 
                      std::numeric_limits<double>::quiet_NaN(), 
                      std::numeric_limits<double>::quiet_NaN()};
    if (verbose)
      printf("i=%d:\tISR(%5.3f,\t%5.3f,\t%.3e,\t%.3e,\t%.3e)\n", i, data_interp[i][0],
              data_interp[i][1], data_interp[i][2], data_interp[i][3], data_interp[i][4]);
  }

  string sname = "_p" + std::to_string(shepard_parameter);
  write_interpolated_time_selected_vals_vtk("SEL" + sname, num_scalar_data,
                                            timestep, data_interp);
  std::fill(data_interp.begin(), data_interp.end(), std::nan);
}


