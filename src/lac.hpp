// 2022-05-18, Juan Carlos Araújo, ju4nk4@gmail.com

#include <stdio.h>
#include <vector>
#include "glm/glm.hpp"
#include "glm/vec3.hpp"
#include "glm/ext.hpp" // printing matrices
#include <glm/gtx/string_cast.hpp>


double sign(double x);
double Lshape (int i, std::vector<double> &x, double xp);
glm::dvec3 rotateAxis (glm::dvec3 R, double angle, glm::dvec3 axis);
glm::dmat3 MatRotateAxis ( double angle, glm::dvec3 axis);
glm::dvec3 rotateAnglesXZ(glm::dvec3 X, double a1, double a2);
glm::dvec3 rotateAnglesZXZ(glm::dvec3 X, double a1, double a2, double a3);
glm::dmat3 matMul (glm::dmat3 a, glm::dmat3 b);

