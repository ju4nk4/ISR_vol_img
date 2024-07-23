// 2022-05-18, Juan Carlos Araújo, ju4nk4@gmail.com

#include "glm/glm.hpp"
#include "glm/vec3.hpp"
#include "glm/ext.hpp" // printing matrices
#include <glm/gtx/string_cast.hpp>

#define GLM_ENABLE_EXPERIMENTAL
#define normalizVec3(val) ((val + glm::dvec3(glm::abs(box.left), glm::abs(box.bottom), glm::abs(box.back))) / normalizedEdges) * glm::dvec3(width, height, depth)

// https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
double sign(double x) {
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

//----------------------------------------------------------------------------------------
// Lagrange interpolants
double Lshape (int i, vector<double> &x, double xp) {
	
	int n = x.size();
	double yp = 1.0;
	
	for(int j=0;j<n;j++) {
		if(i!=j) {
			yp = yp* (xp - x[j])/(x[i] - x[j]);
		}
	}	
	
	return yp;
}


// Rotation with respect to a given axis
// https://followtutorials.com/2012/03/3d-rotation-algorithm-about-arbitrary-axis-with-c-c-code.html
glm::dvec3 rotateAxis (glm::dvec3 R, double angle, glm::dvec3 axis) 
{
    
    double  u = axis.x,
            v = axis.y,
            w = axis.z;
    
    double  L = (u*u + v * v + w * w), sL = glm::sqrt(L),
            u2 = u * u,
            v2 = v * v,
            w2 = w * w,
            sin_an = glm::sin(angle), cos_an = glm::cos(angle);
    
    glm::dmat3 rotationMatrix = glm::dmat3(0.0000);
    
    rotationMatrix[0][0] = (u2 + (v2 + w2) * cos_an) / L;
    rotationMatrix[0][1] = (u * v * (1 - cos_an) - w * sL * sin_an) / L;
    rotationMatrix[0][2] = (u * w * (1 - cos_an) + v * sL * sin_an) / L;
 
    rotationMatrix[1][0] = (u * v * (1 - cos_an) + w * sL * sin_an) / L;
    rotationMatrix[1][1] = (v2 + (u2 + w2) * cos_an) / L;
    rotationMatrix[1][2] = (v * w * (1 - cos_an) - u * sL * sin_an) / L;
 
    rotationMatrix[2][0] = (u * w * (1 - cos_an) - v * sL * sin_an) / L;
    rotationMatrix[2][1] = (v * w * (1 - cos_an) + u * sL * sin_an) / L;
    rotationMatrix[2][2] = (w2 + (u2 + v2) * cos_an) / L;
    
    glm::dvec3 Rnew = rotationMatrix * R;
    
    return Rnew;
}

// returns the matrix
glm::dmat3 MatRotateAxis ( double angle, glm::dvec3 axis) {
    
    double  u = axis.x,
            v = axis.y,
            w = axis.z;
    
    double  L = (u*u + v * v + w * w), sL = glm::sqrt(L),
            u2 = u * u,
            v2 = v * v,
            w2 = w * w,
            sin_an = glm::sin(angle), cos_an = glm::cos(angle);
    
    glm::dmat3 rotationMatrix = glm::dmat3(0.0000);
    
    rotationMatrix[0][0] = (u2 + (v2 + w2) * cos_an) / L;
    rotationMatrix[0][1] = (u * v * (1 - cos_an) - w * sL * sin_an) / L;
    rotationMatrix[0][2] = (u * w * (1 - cos_an) + v * sL * sin_an) / L;
 
    rotationMatrix[1][0] = (u * v * (1 - cos_an) + w * sL * sin_an) / L;
    rotationMatrix[1][1] = (v2 + (u2 + w2) * cos_an) / L;
    rotationMatrix[1][2] = (v * w * (1 - cos_an) - u * sL * sin_an) / L;
 
    rotationMatrix[2][0] = (u * w * (1 - cos_an) - v * sL * sin_an) / L;
    rotationMatrix[2][1] = (v * w * (1 - cos_an) + u * sL * sin_an) / L;
    rotationMatrix[2][2] = (w2 + (u2 + v2) * cos_an) / L;
    
    
    return rotationMatrix;
}


/*

function R = MatRotateAxis(t, mu)
    u = mu(1,1);
    v = mu(2,1);
    w = mu(3,1);
    
    L = (u*u + v * v + w * w);
    u2 = u * u;
    v2 = v * v;
    w2 = w * w;
 
    R(1,1) = (u2 + (v2 + w2) * cos(t)) / L;
    R(1,2) = (u * v * (1 - cos(t)) - w * sqrt(L) * sin(t)) / L;
    R(1,3) = (u * w * (1 - cos(t)) + v * sqrt(L) * sin(t)) / L;
 
    R(2,1) = (u * v * (1 - cos(t)) + w * sqrt(L) * sin(t)) / L;
    R(2,2) = (v2 + (u2 + w2) * cos(t)) / L;
    R(2,3) = (v * w * (1 - cos(t)) - u * sqrt(L) * sin(t)) / L;
 
    R(3,1) = (u * w * (1 - cos(t)) - v * sqrt(L) * sin(t)) / L;
    R(3,2) = (v * w * (1 - cos(t)) + u * sqrt(L) * sin(t)) / L;
    R(3,3) = (w2 + (u2 + v2) * cos(t)) / L;
end

*/

glm::dvec3 rotateAnglesXZ(glm::dvec3 X, double a1, double a2) {
 	
 	glm::dvec3 Y, Y1;
 	glm::dmat3 R = glm::dmat3(0.0000);
 	
 	double t = a1, 	ct = cos(t), 	st = sin(t);
 	
    R[0][0] = 1;
    R[0][1] = 0;
    R[0][2] = 0;
 
    R[1][0] = 0;
    R[1][1] = ct;
    R[1][2] =-st;
 
    R[2][0] = 0;
    R[2][1] = st;
    R[2][2] = ct;
    
    Y1 = R*X;
    
    t = a2; 	ct = cos(t); 	st = sin(t);
 	
 	R[0][0] = ct;
    R[0][1] =-st;
    R[0][2] = 0;
 
    R[1][0] = st;
    R[1][1] = ct;
    R[1][2] = 0;
 
    R[2][0] = 0;
    R[2][1] = 0;
    R[2][2] = 1;
    
    Y = R*Y1;
    
    return Y;
}


glm::dvec3 rotateAnglesZXZ(glm::dvec3 X, double a1, double a2, double a3) {
 	
 	glm::dvec3 Y, Y1, Y2;
 	glm::dmat3 R = glm::dmat3(0.0000);
 	
 	double t = a1, 	ct = cos(t), 	st = sin(t);
 	
 	R[0][0] = ct;
    R[0][1] =-st;
    R[0][2] = 0;
 
    R[1][0] = st;
    R[1][1] = ct;
    R[1][2] = 0;
 
    R[2][0] = 0;
    R[2][1] = 0;
    R[2][2] = 1;
    
    Y1 = R*X;
    
    t = a2; 	ct = cos(t); 	st = sin(t);
 	
    R[0][0] = 1;
    R[0][1] = 0;
    R[0][2] = 0;
 
    R[1][0] = 0;
    R[1][1] = ct;
    R[1][2] =-st;
 
    R[2][0] = 0;
    R[2][1] = st;
    R[2][2] = ct;
    
    Y2 = R*Y1;
    
    t = a3; 	ct = cos(t); 	st = sin(t);
 	
 	R[0][0] = ct;
    R[0][1] =-st;
    R[0][2] = 0;
 
    R[1][0] = st;
    R[1][1] = ct;
    R[1][2] = 0;
 
    R[2][0] = 0;
    R[2][1] = 0;
    R[2][2] = 1;
    
    Y = R*Y2;
    
    return Y;
}

// The implementation of * in glm gives wrong results ... sad!!!
glm::dmat3 matMul (glm::dmat3 a, glm::dmat3 b) {
	
	glm::dmat3 mult(0);
	
	for(unsigned int i = 0; i < 3; ++i)
        for(unsigned int j = 0; j < 3; ++j)
            for(unsigned int k = 0; k < 3; ++k)
            {
                mult[i][j] += a[i][k] * b[k][j];
            }
    return mult;
}

