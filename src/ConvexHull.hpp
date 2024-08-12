// Implementation of convex hull algorithms 
// using the C++ Standard Library.
// Modified by Juan Araújo, ju4nk4@gmail.com, from the code available at:
// https://github.com/MiguelVieira/ConvexHull2D/blob/master/ConvexHull.cpp

#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <vector>



struct point2d {
	float x;
	float y;

	point2d(float xIn, float yIn) : x(xIn), y(yIn) { } 
};

// The z-value of the cross product of segments 
// (a, b) and (a, c). Positive means c is ccw
// from (a, b), negative cw. Zero means its collinear.

float ccw(const point2d& a, const point2d& b, const point2d& c);

// Returns true if a is lexicographically before b.
bool isLeftOf(const point2d& a, const point2d& b) ;


// The length of segment (a, b).
float len(const point2d& a, const point2d& b) ;

// The unsigned distance of p from segment (a, b).
float dist(const point2d& a, const point2d& b, const point2d& p);

// Returns the index of the farthest point from segment (a, b).
size_t getFarthest(const point2d& a, const point2d& b, const std::vector<point2d>& v) ;


// Recursive call of the quickhull algorithm.
void quickHull(const std::vector<point2d>& v, const point2d& a, const point2d& b, 
			   std::vector<point2d>& hull) ;

// QuickHull algorithm. 
// https://en.wikipedia.org/wiki/QuickHull
std::vector<point2d> quickHull(const std::vector<point2d>& v) ;

std::vector<point2d> getPoints() ;

void print(const std::vector<point2d>& v) ;


// Output indices of sorted vector. The vector is not sorted.
void sort_index (std::vector<double> &in_v, std::vector<int> &out_idx) ;

// r: e_polar[i][0], theta := e_polar[i][1]
// vector<double[2]> 
void compute_quickHull_polar (std::vector<std::array<double,2> > e_polar, std::vector<std::array<double,2> > &out) ;



