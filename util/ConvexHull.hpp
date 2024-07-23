// Implementation of convex hull algorithms 
// using the C++ Standard Library.
// Modified by Juan Araújo, ju4nk4@gmail.com, from the code available at:
// https://github.com/MiguelVieira/ConvexHull2D/blob/master/ConvexHull.cpp

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

struct point2d {
	float x;
	float y;

	point2d(float xIn, float yIn) : x(xIn), y(yIn) { } 
};

// The z-value of the cross product of segments 
// (a, b) and (a, c). Positive means c is ccw
// from (a, b), negative cw. Zero means its collinear.

float ccw(const point2d& a, const point2d& b, const point2d& c) {
	return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

// Returns true if a is lexicographically before b.
bool isLeftOf(const point2d& a, const point2d& b) {
	return (a.x < b.x || (a.x == b.x && a.y < b.y));
}


// The length of segment (a, b).
float len(const point2d& a, const point2d& b) {
	return sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
}

// The unsigned distance of p from segment (a, b).
float dist(const point2d& a, const point2d& b, const point2d& p) {
	return fabs((b.x - a.x) * (a.y - p.y) - (b.y - a.y) * (a.x - p.x)) / len(a, b);
}

// Returns the index of the farthest point from segment (a, b).
size_t getFarthest(const point2d& a, const point2d& b, const vector<point2d>& v) {
	size_t idxMax = 0;
	float distMax = dist(a, b, v[idxMax]);

	for (size_t i = 1; i < v.size(); ++i) {
		float distCurr = dist(a, b, v[i]);
		if (distCurr > distMax) {
			idxMax = i;
			distMax = distCurr;
		}
	}

	return idxMax;
}


// Recursive call of the quickhull algorithm.
void quickHull(const vector<point2d>& v, const point2d& a, const point2d& b, 
			   vector<point2d>& hull) {
	if (v.empty()) {
		return;
	}

	point2d f = v[getFarthest(a, b, v)];

	// Collect points to the left of segment (a, f)
	vector<point2d> left;
	for (auto p : v) {
		if (ccw(a, f, p) > 0) {
			left.push_back(p);
		}
	}
	quickHull(left, a, f, hull);
	
	// Add f to the hull
	hull.push_back(f);

	// Collect points to the left of segment (f, b)
	vector<point2d> right;
	for (auto p : v) {
		if (ccw(f, b, p) > 0) {
			right.push_back(p);
		}
	}
	quickHull(right, f, b, hull);
}

// QuickHull algorithm. 
// https://en.wikipedia.org/wiki/QuickHull
vector<point2d> quickHull(const vector<point2d>& v) {
	vector<point2d> hull;
	
	// Start with the leftmost and rightmost points.
	point2d a = *min_element(v.begin(), v.end(), isLeftOf);
	point2d b = *max_element(v.begin(), v.end(), isLeftOf);

	// Split the points on either side of segment (a, b)
	vector<point2d> left, right;
	for (auto p : v) {
		ccw(a, b, p) > 0 ? left.push_back(p) : right.push_back(p);
	}
	
	// Be careful to add points to the hull
	// in the correct order. Add our leftmost point.
	hull.push_back(a);

	// Add hull points from the left (top)
	quickHull(left, a, b, hull);

	// Add our rightmost point
	hull.push_back(b);

	// Add hull point2ds from the right (bottom)
	quickHull(right, b, a, hull);

	return hull;
}

vector<point2d> getPoints() {
	vector<point2d> v;
	
	const float lo = -100.0;
	const float hi = 100.0;

	for (int i = 0; i < 100; ++i) {
		float x = lo + 
			static_cast<float>(
				rand()) / static_cast<float>(RAND_MAX / (hi - lo));

		float y = lo + 
			static_cast<float>(
				rand()) / static_cast<float>(RAND_MAX / (hi - lo));

		v.push_back(point2d(x, y));
	}

	return v;
}

void print(const vector<point2d>& v) {
	for (auto p : v) {
		cout << "\t" << p.x << ", " << p.y << endl;
	}
}


// Output indices of sorted vector. The vector is not sorted.
void sort_index (vector<double> &in_v, vector<int> &out_idx) {
	//vector<int> V(in_v.size() );
	// https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
	
	out_idx.resize( in_v.size() );
	std::iota(out_idx.begin(),out_idx.end(),0); //Initializing
	sort( out_idx.begin(),out_idx.end(), [&](int i,int j){return in_v[i]<in_v[j];} );
}

// r: e_polar[i][0], theta := e_polar[i][1]
// vector<double[2]> 
void compute_quickHull_polar (vector<std::array<double,2> > e_polar, vector<std::array<double,2> > &out) {
	
	unsigned int e_N = e_polar.size(), N;
	vector<point2d> v;
	
	for (unsigned int i = 0; i < e_N; i++) {
		double x = e_polar[i][0]*cos(e_polar[i][1]), y = e_polar[i][0]*sin(e_polar[i][1]);
			v.push_back( point2d(x,y) );
	}
	vector<point2d> h = quickHull(v);
	
	N = h.size();
	
	vector<double> rv(N), ov(N);
	for (unsigned int i = 0; i < N; i++) {
		double x = (double)(h[i].x), y = (double)(h[i].y), r = sqrt (x*x + y*y), o = atan2(y,x);
		rv[i] = r;
		ov[i] = o;
	}
	
	vector<int> idx;
	sort_index (ov, idx);
	
	out.resize (N);
	for (unsigned int i = 0; i < N; i++) {
		out[i][0] = rv[ idx[i] ];
		out[i][1] = ov[ idx[i] ];
	}
	
	//return out;
}




