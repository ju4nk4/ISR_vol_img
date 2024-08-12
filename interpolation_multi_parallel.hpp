#include <array>
//#include <glm/glm.hpp>
//#include <string>
#include <vector>

void interpMultiVecParallel(float p, unsigned nd, const unsigned n_scalars,
                            const float* xd, const float* fd, unsigned ni,
                            float* xi, bool estimate_derivatives, float c_nu,
                            const unsigned NC, float* Fi);

void spInterpMultiVecParallel(const unsigned ND, const unsigned NS,
                              const unsigned nx, const unsigned ny,
                              const unsigned nz, const unsigned NI, float* xi,
                              const float* const data_array, const unsigned NC,
                              float* Fi);

float* allocateAlignedResult(const unsigned NI, const unsigned NS,
                             const unsigned VEC_WIDTH);

void unpackInterpolation(const unsigned NI, const unsigned NS,
                         const unsigned VEC_WIDTH, float* Fi,
                         std::vector<std::array<double, 5>>& result);
