#include "interpolation_multi_parallel.hpp"

#include <immintrin.h>

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <glm/glm.hpp>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "omp.h"

/*
 * Terminology
 * IP: Interpolation Point
 * DP: Data Point
 */

void derivativesMultiVecParallel(const unsigned NS, const unsigned dim,
                                 const unsigned ND, const unsigned PADDED_ND,
                                 const float c_nu, const float* const xd,
                                 const float* const fd, const unsigned NC,
                                 float* Fxyz, float* min_f, float* max_f,
                                 float* nu) {
  const unsigned VEC_WIDTH = 8U;
  const unsigned SCALAR_STRIDE = VEC_WIDTH;
  const unsigned BLOCK_STRIDE = SCALAR_STRIDE * NS;
  const unsigned NUM_WEIGHTS = VEC_WIDTH * NS;
  const unsigned AXIS_STRIDE = PADDED_ND * NS;

  for (unsigned k = 0; k < NS; k++) {
    min_f[k] = 1.0e20f;
    max_f[k] = -1.0e20f;

    // find min and max
    for (unsigned i = 0; i < ND; i++) {
      if (!std::isnan(fd[ND * k + i])) {
        min_f[k] = std::min(min_f[k], fd[ND * k + i]);
        max_f[k] = std::max(max_f[k], fd[ND * k + i]);
      }
    }
  }

  const unsigned NUM_BLOCKS = ND / VEC_WIDTH;
  const unsigned LEFTOVERS = ND % VEC_WIDTH;

  omp_set_num_threads(NC);
  float max_Fx2[NS * NC];
  float max_Fy2[NS * NC];
  float max_Fz2[NS * NC];

  for (unsigned i = 0; i < NS * NC; i++) {
    max_Fx2[i] = 0.0f;
    max_Fy2[i] = 0.0f;
    max_Fz2[i] = 0.0f;
  }

#pragma omp parallel default(shared)
  {
    unsigned thread_id = omp_get_thread_num();
    float S[NUM_WEIGHTS];
    const float* base_xd_i;
    const float* base_xd_j;
    const float* base_fd_i;
    const float* base_fd_j;
    float* base_Fxyz;

    float max_Fx2_local[NS];
    float max_Fy2_local[NS];
    float max_Fz2_local[NS];

    for (unsigned i = 0; i < NS; i++) {
      max_Fx2_local[i] = 0.0f;
      max_Fy2_local[i] = 0.0f;
      max_Fz2_local[i] = 0.0f;
    }

    float dx, dy, dz, t2, wji, WxdF_t2;

    /****************** VECTOR REGISTERS ******************/
    __m256 v_x_i;    // x-coordinate block-i points
    __m256 v_y_i;    // y-coordinate block-i points
    __m256 v_z_i;    // z-coordinate block-i points
    __m256 v_x_j;    // x-coordinate block-j points
    __m256 v_y_j;    // y-coordinate block-j points
    __m256 v_z_j;    // z-coordinate block-j points
    __m256 v_dx;     // x-coordinate difference
    __m256 v_dy;     // y-coordinate difference
    __m256 v_dz;     // z-coordinate difference
    __m256 v_dji;    // distance between each pair of points
    __m256 v_fd_i;   // scalar values of block-i points
    __m256 v_fd_j;   // scalar values of block_j points
    __m256 v_Wji;    // weights of each pair of points
    __m256 v_Wji_m;  // weights after mask
    __m256 v_Wji2;   // square of the weights
    __m256 v_S;      // Sum of weights for each point
    __m256 v_mask;   // combined mask for i and j
    __m256 v_Fx;     //
    __m256 v_Fy;     //
    __m256 v_Fz;     //

    __m256 ones = _mm256_set1_ps(1.0f);
    __m256i cycle_lower = _mm256_set_epi32(
        0, 7, 6, 5, 4, 3, 2, 1);  // Cycle-shift every float to a lower
                                  // address, with lowest going into highest
    __m256i pattern = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

#pragma omp for schedule(static)
    for (unsigned BLOCK_I = 0; BLOCK_I < NUM_BLOCKS; BLOCK_I++) {
      base_xd_i = xd + BLOCK_I * VEC_WIDTH;
      base_fd_i = fd + BLOCK_I * VEC_WIDTH;
      base_Fxyz = Fxyz + BLOCK_I * BLOCK_STRIDE;

      for (unsigned i = 0; i < NUM_WEIGHTS; i++) {
        S[i] = 0.0f;
        base_Fxyz[i] = 0.0f;
        base_Fxyz[AXIS_STRIDE + i] = 0.0f;
        base_Fxyz[2U * AXIS_STRIDE + i] = 0.0f;
      }

      v_x_i = _mm256_loadu_ps(base_xd_i);
      v_y_i = _mm256_loadu_ps(base_xd_i + ND);
      v_z_i = _mm256_loadu_ps(base_xd_i + 2U * ND);

      for (unsigned BLOCK_J = 0; BLOCK_J < NUM_BLOCKS; BLOCK_J++) {
        base_xd_j = xd + BLOCK_J * VEC_WIDTH;
        base_fd_j = fd + BLOCK_J * VEC_WIDTH;

        v_x_j = _mm256_loadu_ps(base_xd_j);
        v_y_j = _mm256_loadu_ps(base_xd_j + ND);
        v_z_j = _mm256_loadu_ps(base_xd_j + 2U * ND);

        // modify this variable if BLOCK_I == BLOCK_J to avoid first computation
        // and cycle-shift the corresponding values.
        unsigned REUSE_ITER = 0U;

        if (BLOCK_I == BLOCK_J) {
          v_x_j = _mm256_permutevar8x32_ps(v_x_j, cycle_lower);
          v_y_j = _mm256_permutevar8x32_ps(v_y_j, cycle_lower);
          v_z_j = _mm256_permutevar8x32_ps(v_z_j, cycle_lower);
          pattern = _mm256_permutevar8x32_epi32(pattern, cycle_lower);
          REUSE_ITER++;
        }

        // FOR LOOP TO REUSE DATA
        for (; REUSE_ITER < VEC_WIDTH; REUSE_ITER++) {
          v_dx = _mm256_sub_ps(v_x_j, v_x_i);
          v_dji = _mm256_mul_ps(v_dx, v_dx);

          v_dy = _mm256_sub_ps(v_y_j, v_y_i);
          v_dji = _mm256_fmadd_ps(v_dy, v_dy, v_dji);

          v_dz = _mm256_sub_ps(v_z_j, v_z_i);
          v_dji = _mm256_fmadd_ps(v_dz, v_dz, v_dji);

          // Cycle-shift data points' coordinates
          v_x_j = _mm256_permutevar8x32_ps(v_x_j, cycle_lower);
          v_y_j = _mm256_permutevar8x32_ps(v_y_j, cycle_lower);
          v_z_j = _mm256_permutevar8x32_ps(v_z_j, cycle_lower);

          v_Wji = _mm256_div_ps(ones, v_dji);
          v_Wji2 = _mm256_mul_ps(v_Wji, v_Wji);

          // FOR LOOP FOR THE DIFFERENT SCALARS
          for (unsigned k = 0; k < NS; k++) {
            v_fd_i = _mm256_loadu_ps(base_fd_i + ND * k);
            v_fd_j = _mm256_loadu_ps(base_fd_j + ND * k);
            v_fd_j = _mm256_permutevar8x32_ps(v_fd_j, pattern);

            v_mask = _mm256_cmp_ps(v_fd_i, v_fd_j, _CMP_ORD_Q);
            v_fd_i = _mm256_and_ps(v_fd_i, v_mask);
            v_fd_j = _mm256_and_ps(v_fd_j, v_mask);

            v_fd_i = _mm256_sub_ps(v_fd_j, v_fd_i);
            v_fd_i = _mm256_mul_ps(v_fd_i, v_Wji2);

            v_Wji_m = _mm256_and_ps(v_Wji, v_mask);
            v_S = _mm256_loadu_ps(S + VEC_WIDTH * k);
            v_S = _mm256_add_ps(v_S, v_Wji_m);
            _mm256_storeu_ps(S + VEC_WIDTH * k, v_S);

            v_Fx = _mm256_load_ps(base_Fxyz + SCALAR_STRIDE * k);
            v_Fx = _mm256_fmadd_ps(v_fd_i, v_dx, v_Fx);
            _mm256_store_ps(base_Fxyz + SCALAR_STRIDE * k, v_Fx);

            v_Fy = _mm256_load_ps(base_Fxyz + AXIS_STRIDE + SCALAR_STRIDE * k);
            v_Fy = _mm256_fmadd_ps(v_fd_i, v_dy, v_Fy);
            _mm256_store_ps(base_Fxyz + AXIS_STRIDE + SCALAR_STRIDE * k, v_Fy);

            v_Fz = _mm256_load_ps(base_Fxyz + AXIS_STRIDE * 2U +
                                  SCALAR_STRIDE * k);
            v_Fz = _mm256_fmadd_ps(v_fd_i, v_dz, v_Fz);
            _mm256_store_ps(base_Fxyz + AXIS_STRIDE * 2U + SCALAR_STRIDE * k,
                            v_Fz);
          }  // LOOP SCALARS

          pattern = _mm256_permutevar8x32_epi32(pattern, cycle_lower);
        }  // LOOP REUSE DATA
      }  // LOOP BLOCK_J

      // HANDLE LEFTOVERS HERE
      // @todo: vectorize this
      base_xd_j = xd + NUM_BLOCKS * VEC_WIDTH;
      base_fd_j = fd + NUM_BLOCKS * VEC_WIDTH;
      for (unsigned i = 0; i < VEC_WIDTH; i++) {
        for (unsigned j = 0; j < LEFTOVERS; j++) {
          dx = base_xd_j[j] - base_xd_i[i];
          dy = base_xd_j[j + ND] - base_xd_i[i + ND];
          dz = base_xd_j[j + 2U * ND] - base_xd_i[i + 2U * ND];
          t2 = dx * dx + dy * dy + dz * dz;

          wji = 1.0f / t2;
          for (unsigned k = 0; k < NS; k++) {
            if (!std::isnan(base_fd_j[j + ND * k]) and
                !std::isnan(base_fd_i[i + ND * k])) {
              S[SCALAR_STRIDE * k + i] += wji;
              WxdF_t2 =
                  wji * wji * (base_fd_j[j + ND * k] - base_fd_i[i + ND * k]);
              base_Fxyz[SCALAR_STRIDE * k + i] += WxdF_t2 * dx;
              base_Fxyz[AXIS_STRIDE + SCALAR_STRIDE * k + i] += WxdF_t2 * dy;
              base_Fxyz[AXIS_STRIDE * 2U + SCALAR_STRIDE * k + i] +=
                  WxdF_t2 * dz;
            }
          }
        }
      }

      for (unsigned k = 0; k < NS; k++) {
        for (unsigned ii = 0; ii < VEC_WIDTH; ii++) {
          if (S[ii + VEC_WIDTH * k] > 1e-7f) {
            base_Fxyz[SCALAR_STRIDE * k + ii] /= S[ii + VEC_WIDTH * k];
            base_Fxyz[AXIS_STRIDE + SCALAR_STRIDE * k + ii] /=
                S[ii + VEC_WIDTH * k];
            base_Fxyz[AXIS_STRIDE * 2U + SCALAR_STRIDE * k + ii] /=
                S[ii + VEC_WIDTH * k];
          }

          max_Fx2_local[k] =
              std::max(max_Fx2_local[k], base_Fxyz[SCALAR_STRIDE * k + ii] *
                                             base_Fxyz[SCALAR_STRIDE * k + ii]);
          max_Fy2_local[k] =
              std::max(max_Fy2_local[k],
                       base_Fxyz[AXIS_STRIDE + SCALAR_STRIDE * k + ii] *
                           base_Fxyz[AXIS_STRIDE + SCALAR_STRIDE * k + ii]);
          max_Fz2_local[k] = std::max(
              max_Fz2_local[k],
              base_Fxyz[AXIS_STRIDE * 2U + SCALAR_STRIDE * k + ii] *
                  base_Fxyz[AXIS_STRIDE * 2U + SCALAR_STRIDE * k + ii]);
        }
      }  // LOOP SCALARS
    }  // LOOP BLOCK_I

    for (unsigned k = 0; k < NS; k++) {
      max_Fx2[k + thread_id * NS] = max_Fx2_local[k];
      max_Fy2[k + thread_id * NS] = max_Fy2_local[k];
      max_Fz2[k + thread_id * NS] = max_Fz2_local[k];
    }
  }  // END PARALLEL REGION

  // reduction of max -- only for master thread
  for (unsigned k = 0; k < NS; k++) {
    for (unsigned th_id = 1; th_id < NC; th_id++) {
      max_Fx2[k] = std::max(max_Fx2[k], max_Fx2[k + NS * th_id]);
      max_Fy2[k] = std::max(max_Fy2[k], max_Fy2[k + NS * th_id]);
      max_Fz2[k] = std::max(max_Fz2[k], max_Fz2[k + NS * th_id]);
    }
  }

  // LEFTOVERS END OF I
  float* base_Fxyz = Fxyz + BLOCK_STRIDE * NUM_BLOCKS;
  const float* base_xd_i = xd + NUM_BLOCKS * VEC_WIDTH;
  const float* base_fd_i = fd + NUM_BLOCKS * VEC_WIDTH;
  float dx, dy, dz, t2, wji, WxdF_t2;
  float s[NS];
  for (unsigned i = 0; i < LEFTOVERS; i++) {
    for (unsigned k = 0; k < NS; k++) {
      s[k] = 0.0f;
      base_Fxyz[i + VEC_WIDTH * k] = 0.0f;
      base_Fxyz[i + VEC_WIDTH * k + AXIS_STRIDE] = 0.0f;
      base_Fxyz[i + VEC_WIDTH * k + AXIS_STRIDE * 2U] = 0.0f;
    }
    for (unsigned j = 0; j < ND; j++) {
      if (j != (i + NUM_BLOCKS * VEC_WIDTH)) {
        dx = xd[j] - base_xd_i[i];
        dy = xd[j + ND] - base_xd_i[i + ND];
        dz = xd[j + ND * 2U] - base_xd_i[i + ND * 2U];
        t2 = dx * dx + dy * dy + dz * dz;
        wji = 1.0f / t2;
        for (unsigned k = 0; k < NS; k++) {
          if (!std::isnan(fd[j + ND * k]) and
              !std::isnan(base_fd_i[i + ND * k])) {
            s[k] += wji;
            WxdF_t2 = wji * wji * (fd[j + ND * k] - base_fd_i[i + ND * k]);
            base_Fxyz[i + SCALAR_STRIDE * k] += WxdF_t2 * dx;
            base_Fxyz[i + SCALAR_STRIDE * k + AXIS_STRIDE] += WxdF_t2 * dy;
            base_Fxyz[i + SCALAR_STRIDE * k + AXIS_STRIDE * 2U] += WxdF_t2 * dz;
          }
        }
      }
    }
    for (unsigned k = 0; k < NS; k++) {
      if (s[k] > 1e-7) {
        base_Fxyz[i + SCALAR_STRIDE * k] /= s[k];
        base_Fxyz[i + SCALAR_STRIDE * k + AXIS_STRIDE] /= s[k];
        base_Fxyz[i + SCALAR_STRIDE * k + AXIS_STRIDE * 2U] /= s[k];
      }

      max_Fx2[k] = std::max(max_Fx2[k], base_Fxyz[i + SCALAR_STRIDE * k] *
                                            base_Fxyz[i + SCALAR_STRIDE * k]);
      max_Fy2[k] = std::max(max_Fy2[k],
                            base_Fxyz[i + SCALAR_STRIDE * k + AXIS_STRIDE] *
                                base_Fxyz[i + SCALAR_STRIDE * k + AXIS_STRIDE]);
      max_Fz2[k] = std::max(
          max_Fz2[k], base_Fxyz[i + SCALAR_STRIDE * k + AXIS_STRIDE * 2U] *
                          base_Fxyz[i + SCALAR_STRIDE * k + AXIS_STRIDE * 2U]);
    }
  }

  for (unsigned k = 0; k < NS; k++) {
    if (max_Fx2[k] + max_Fy2[k] + max_Fz2[k] > 0.0f) {
      nu[k] = c_nu * (max_f[k] - min_f[k]) /
              sqrt(max_Fx2[k] + max_Fy2[k] + max_Fz2[k]);
    } else {
      nu[k] = 0.0f;
    }

    // std::cout << "max_Fx2[" << k << "] = " << max_Fx2[k] << std::endl;
    // std::cout << "max_Fy2[" << k << "] = " << max_Fy2[k] << std::endl;
    // std::cout << "max_Fz2[" << k << "] = " << max_Fz2[k] << std::endl;
  }
}

/**
 * @brief
 *
 * @param p            is there any difference between p and dim?
 * @param ND           number of input DPs
 * @param NS           number of independent scalars (5 in this case)
 * @param xd           array with coordinates of DP
 * @param fd           scalars of DP
 * @param NI           number of IP
 * @param xi           array with the coordinates of the IP
 * @param estimate_derivatives
 * @param c_nu         double
 * @param NC           number of cores used in the computation
 * @param Fi           Interpolated scalars (the result)
 */
void interpMultiVecParallel(float p, const unsigned ND, const unsigned NS,
                            const float* xd, const float* fd, unsigned NI,
                            float* xi, bool estimate_derivatives, float c_nu,
                            const unsigned NC, float* Fi) {
  constexpr unsigned dim = 3U;
  // constexpr unsigned sf = 3;  // @warning: what is sf?

  unsigned VEC_WIDTH = 8U;
  const unsigned SCALAR_STRIDE = VEC_WIDTH;
  const unsigned BLOCK_STRIDE = NS * VEC_WIDTH;  // for non-consecutive xyz
  const unsigned ELEMS_PER_BLOCK = NS * VEC_WIDTH;

  const unsigned PADDED_ND =
      (ND % VEC_WIDTH == 0) ? ND : ((ND / VEC_WIDTH) + 1) * VEC_WIDTH;
  const unsigned AXIS_STRIDE = PADDED_ND * NS;

  const unsigned alignment = VEC_WIDTH * sizeof(float);
  const unsigned size_Fxyz = 3U * PADDED_ND * NS;
  float* Fxyz = (float*)aligned_alloc(alignment, size_Fxyz * sizeof(float));

  float min_f[NS];
  float max_f[NS];
  float nu[NS];

  // auto pre_der = std::chrono::high_resolution_clock::now();
  derivativesMultiVecParallel(NS, dim, ND, PADDED_ND, c_nu, xd, fd, NC, Fxyz,
                              min_f, max_f, nu);

  // --------------------------------------------------------------------------
  // ---------------------------- BEGIN PARALLELIZE ---------------------------
  // --------------------------------------------------------------------------

  unsigned N_BLOCK_INTERP = NI / VEC_WIDTH;
  unsigned N_BLOCK_DATA = ND / VEC_WIDTH;
  unsigned LEFTOVERS = ND % VEC_WIDTH;
  omp_set_num_threads(NC);

#pragma omp parallel default(shared)
  {
    float S[ELEMS_PER_BLOCK];
    float Fi_local[ELEMS_PER_BLOCK];
    const float* base_xi;
    const float* base_xd;
    const float* base_fd;
    const float* base_Fxyz;
    float* base_Fi;

    /****************** VECTOR REGISTERS ******************/
    __m256 v_xi_x;   // x-coordinate interpolation points
    __m256 v_xi_y;   // y-coordinate interpolation points
    __m256 v_xi_z;   // z-coordinate interpolation points
    __m256 v_xd_x;   // x-coordinate data points
    __m256 v_xd_y;   // y-coordinate data points
    __m256 v_xd_z;   // z-coordinate data points
    __m256 v_dx;     // x-coordinate difference
    __m256 v_dy;     // y-coordinate difference
    __m256 v_dz;     // z-coordinate difference
    __m256 v_dij;    // distance between each pair of points
    __m256 v_dij_r;  // reciprocal of v_dij
    __m256 v_W;      // weights of each pair of points
    __m256 v_W_m;    // weights after mask
    __m256 v_mask;   // mask for nans
    __m256 v_S;      // Sum of weights for each interpolation point
    __m256 v_Fi;     // Scalar values in interpolation points
    __m256 v_fd;     // Scalar values in data points
    __m256 v_FX;     // x-coordinate scalar derivative in data points
    __m256 v_FY;     // y-coordinate scalar derivative in data points
    __m256 v_FZ;     // z-coordinate scalar derivative in data points
    __m256 v_nu;     // value of nu for the corresponding scalar
    __m256 ones = _mm256_set1_ps(1.0f);
    __m256i cycle_lower = _mm256_set_epi32(
        0, 7, 6, 5, 4, 3, 2, 1);  // Cycle-shift every float to a lower
                                  // address, with lowest going into highest

// BI_ID = Block-Interpolation Identifier
#pragma omp for schedule(static)
    for (unsigned BI_ID = 0; BI_ID < N_BLOCK_INTERP; BI_ID++) {
      // point to the correct point of the x coordinate.
      base_xi = xi + VEC_WIDTH * BI_ID;
      base_Fi = Fi + ELEMS_PER_BLOCK * BI_ID;

      // Loop to initialize data
      for (unsigned i = 0; i < ELEMS_PER_BLOCK; i++) {
        S[i] = 0.0f;
        Fi_local[i] = 0.0f;
      }

      v_xi_x = _mm256_loadu_ps(base_xi);
      v_xi_y = _mm256_loadu_ps(base_xi + NI);
      v_xi_z = _mm256_loadu_ps(base_xi + 2U * NI);

      for (unsigned BD_ID = 0; BD_ID < N_BLOCK_DATA; BD_ID++) {
        base_xd = xd + VEC_WIDTH * BD_ID;
        base_fd = fd + VEC_WIDTH * BD_ID;
        base_Fxyz = Fxyz + BLOCK_STRIDE * BD_ID;

        v_xd_x = _mm256_loadu_ps(base_xd);
        v_xd_y = _mm256_loadu_ps(base_xd + ND);
        v_xd_z = _mm256_loadu_ps(base_xd + 2U * ND);

        // FOR LOOP TO REUSE DATA
        for (unsigned j_BLOCK = 0; j_BLOCK < VEC_WIDTH; j_BLOCK++) {
          v_dx = _mm256_sub_ps(v_xi_x, v_xd_x);
          v_dy = _mm256_sub_ps(v_xi_y, v_xd_y);
          v_dz = _mm256_sub_ps(v_xi_z, v_xd_z);

          // Cycle-shift data points' coordinates for next iteration
          v_xi_x = _mm256_permutevar8x32_ps(v_xi_x, cycle_lower);
          v_xi_y = _mm256_permutevar8x32_ps(v_xi_y, cycle_lower);
          v_xi_z = _mm256_permutevar8x32_ps(v_xi_z, cycle_lower);

          v_dij = _mm256_mul_ps(v_dx, v_dx);
          v_dij = _mm256_fmadd_ps(v_dy, v_dy, v_dij);
          v_dij = _mm256_fmadd_ps(v_dz, v_dz, v_dij);

          v_dij = _mm256_sqrt_ps(v_dij);
          v_dij_r = _mm256_div_ps(ones, v_dij);
          v_W = _mm256_mul_ps(v_dij_r, v_dij_r);
          v_W = _mm256_mul_ps(v_W, v_dij_r);

          // FOR LOOP FOR THE DIFFERENT SCALARS
          for (unsigned k = 0; k < NS; k++) {
            v_fd = _mm256_loadu_ps(base_fd + ND * k);

            v_mask = _mm256_cmp_ps(v_fd, v_fd, _CMP_ORD_Q);
            v_fd = _mm256_and_ps(v_fd, v_mask);
            v_W_m = _mm256_and_ps(v_W, v_mask);

            v_S = _mm256_loadu_ps(S + VEC_WIDTH * k);
            v_S = _mm256_add_ps(v_S, v_W_m);
            v_S = _mm256_permutevar8x32_ps(v_S, cycle_lower);
            _mm256_storeu_ps(S + VEC_WIDTH * k, v_S);

            // Use v_FX as accumulator
            v_FX = _mm256_load_ps(base_Fxyz + SCALAR_STRIDE * k);
            v_FX = _mm256_mul_ps(v_FX, v_dx);

            v_FY = _mm256_load_ps(base_Fxyz + AXIS_STRIDE + SCALAR_STRIDE * k);
            v_FX = _mm256_fmadd_ps(v_FY, v_dy, v_FX);

            v_FZ = _mm256_load_ps(base_Fxyz + 2U * AXIS_STRIDE +
                                  SCALAR_STRIDE * k);
            v_FX = _mm256_fmadd_ps(v_FZ, v_dz, v_FX);

            v_nu = _mm256_broadcast_ss(nu + k);
            v_nu = _mm256_div_ps(v_nu, _mm256_add_ps(v_nu, v_dij));

            v_FX = _mm256_mul_ps(v_FX, v_nu);

            v_Fi = _mm256_loadu_ps(Fi_local + VEC_WIDTH * k);
            v_Fi = _mm256_fmadd_ps(v_W_m, v_fd, v_Fi);
            v_Fi = _mm256_fmadd_ps(v_FX, v_W_m, v_Fi);
            v_Fi = _mm256_permutevar8x32_ps(v_Fi, cycle_lower);
            _mm256_storeu_ps(Fi_local + VEC_WIDTH * k, v_Fi);
          }  // LOOP SCALARS

        }  // LOOP DIFFERENT COMBINATIONS
      }  // LOOP BLOCK DATA POINTS

      // ACCESSED IF ND IS NOT A MULTIPLE OF VEC_WIDTH
      base_xd = xd + VEC_WIDTH * N_BLOCK_DATA;
      base_fd = fd + VEC_WIDTH * N_BLOCK_DATA;
      base_Fxyz = Fxyz + ELEMS_PER_BLOCK * N_BLOCK_DATA;
      for (unsigned j = 0; j < LEFTOVERS; j++) {
        v_xd_x = _mm256_broadcast_ss(base_xd + j);
        v_xd_y = _mm256_broadcast_ss(base_xd + ND + j);
        v_xd_z = _mm256_broadcast_ss(base_xd + ND * 2U + j);

        v_dx = _mm256_sub_ps(v_xi_x, v_xd_x);
        v_dy = _mm256_sub_ps(v_xi_y, v_xd_y);
        v_dz = _mm256_sub_ps(v_xi_z, v_xd_z);

        v_dij = _mm256_mul_ps(v_dx, v_dx);
        v_dij = _mm256_fmadd_ps(v_dy, v_dy, v_dij);
        v_dij = _mm256_fmadd_ps(v_dz, v_dz, v_dij);

        v_dij = _mm256_sqrt_ps(v_dij);
        v_dij_r = _mm256_div_ps(ones, v_dij);
        v_W = _mm256_mul_ps(v_dij_r, v_dij_r);
        v_W = _mm256_mul_ps(v_W, v_dij_r);

        for (unsigned k = 0; k < NS; k++) {
          v_fd = _mm256_broadcast_ss(base_fd + ND * k + j);
          v_mask = _mm256_cmp_ps(v_fd, v_fd, _CMP_ORD_Q);
          v_fd = _mm256_and_ps(v_fd, v_mask);
          v_W_m = _mm256_and_ps(v_W, v_mask);

          v_S = _mm256_loadu_ps(S + VEC_WIDTH * k);
          v_S = _mm256_add_ps(v_S, v_W_m);
          _mm256_storeu_ps(S + VEC_WIDTH * k, v_S);

          v_Fi = _mm256_loadu_ps(Fi_local + VEC_WIDTH * k);
          v_Fi = _mm256_fmadd_ps(v_W_m, v_fd, v_Fi);

          v_FX = _mm256_broadcast_ss(base_Fxyz + SCALAR_STRIDE * k + j);
          v_FX = _mm256_mul_ps(v_FX, v_dx);

          v_FY = _mm256_broadcast_ss(base_Fxyz + AXIS_STRIDE +
                                     SCALAR_STRIDE * k + j);
          v_FX = _mm256_fmadd_ps(v_FY, v_dy, v_FX);

          v_FZ = _mm256_broadcast_ss(base_Fxyz + 2U * AXIS_STRIDE +
                                     SCALAR_STRIDE * k + j);
          v_FX = _mm256_fmadd_ps(v_FZ, v_dz, v_FX);

          v_nu = _mm256_broadcast_ss(nu + k);
          v_nu = _mm256_div_ps(v_nu, _mm256_add_ps(v_nu, v_dij));

          v_FX = _mm256_mul_ps(v_FX, v_nu);
          v_Fi = _mm256_fmadd_ps(v_FX, v_W_m, v_Fi);
          _mm256_storeu_ps(Fi_local + VEC_WIDTH * k, v_Fi);
        }
      }

      for (unsigned k = 0; k < NS; k++) {
        v_Fi = _mm256_loadu_ps(Fi_local + VEC_WIDTH * k);
        v_S = _mm256_load_ps(S + VEC_WIDTH * k);
        v_Fi = _mm256_div_ps(v_Fi, v_S);
        _mm256_store_ps(base_Fi + VEC_WIDTH * k, v_Fi);

        // @todo: add comparison with min_f?
      }

    }  // LOOP BLOCK INTERPOLATION POINTS
  }  // PARALLEL REGION

  // --------------------------------------------------------------------------
  // ----------------------------- END PARALLELIZE ----------------------------
  // --------------------------------------------------------------------------

  free(Fxyz);
}

// Preparing data for volume interpolation for NEL data
void spInterpMultiVecParallel(const unsigned ND, const unsigned NS,
                              const unsigned nx, const unsigned ny,
                              const unsigned nz, const unsigned NI, float* xi,
                              const float* const data_array, const unsigned NC,
                              float* Fi) {
  bool estimate_derivatives = true;

  // @warning: What is p?
  float c_nu = 0.5, p = 3.0;

  float* xd = new float[ND * 3U];
  float* fd = new float[ND * NS];

  // @warning: COPY OF DATA
  for (unsigned i = 0; i < ND; i++) {
    xd[0 * ND + i] = data_array[i * (3 + NS) + 0];
    xd[1 * ND + i] = data_array[i * (3 + NS) + 1];
    xd[2 * ND + i] = data_array[i * (3 + NS) + 2];

    for (unsigned j = 0; j < NS; j++) {
      fd[j * ND + i] = data_array[i * (3 + NS) + 3 + j];
    }
  }

  auto begin = std::chrono::high_resolution_clock::now();
  interpMultiVecParallel(p, ND, NS, xd, fd, NI, xi, estimate_derivatives, c_nu,
                         NC, Fi);
  auto end = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration<double>(end - begin).count();
  std::cout << "Time multi-scalar interpolation: " << elapsed << " seconds.\n";

  delete[] xd;
  delete[] fd;
}

float* allocateAlignedResult(const unsigned NI, const unsigned NS,
                             const unsigned VEC_WIDTH) {
  const unsigned alignment = VEC_WIDTH * sizeof(float);
  const unsigned size_Fi = NI * NS;
  float* Fi = (float*)aligned_alloc(alignment, size_Fi * sizeof(float));
  return Fi;
}

void unpackInterpolation(const unsigned NI, const unsigned NS,
                         const unsigned VEC_WIDTH, float* Fi,
                         std::vector<std::array<double, 5>>& result) {
  unsigned BLOCK_STRIDE = NS * VEC_WIDTH;

  for (unsigned i = 0; i < NI; i++) {
    unsigned blck_id = i / VEC_WIDTH;
    unsigned in_blck_id = i % VEC_WIDTH;
    for (unsigned k = 0; k < NS; k++) {
      result[i][k] = Fi[blck_id * BLOCK_STRIDE + VEC_WIDTH * k + in_blck_id];
    }
  }
}
