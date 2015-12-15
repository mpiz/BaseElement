#pragma once

#include "macro.h"

#include <fstream>
#include <cmath>
#include <utility>


using namespace std;

class node {
 public:
	 node();
	 node(double s_x, double s_y, double s_z, int num = 0);
	 node (double s_n[3]);
	 double& operator [] (int i);
	 node& operator = (const node pr);
	 bool operator < (const node pr);
	 bool operator > (const node pr);
	 friend ostream& operator << (ostream& os, node& out_n);
	 friend istream& operator >> (istream& is, node& inp_n);

	 double x, y, z;
	 int number;
};

class point {
 public:
	 point();
	 point(double s_x, double s_y, double s_z);
	 point (double s_n[3]);
	 point(node nd);
	 double& operator [] (int i);
	 point& operator = (const point pr);
	 point& operator = (const node pr);

	 double x, y, z;
};

class vec3d {

 public:
	vec3d();
	vec3d(double s_x, double s_y, double s_z);
	vec3d(point start, point end);
	vec3d(node start, node end);
	vec3d(node nd);

	double& operator [] (int i);
	double operator * (const vec3d pr);
	vec3d operator + (const vec3d pr);
	vec3d operator - (const vec3d pr);
	bool operator == (const vec3d pr);
	bool operator != (const vec3d pr);
	friend vec3d operator * (const double a, const vec3d vec);
	vec3d operator / (const double a);
	friend vec3d operator * (const matrix(3) M, vec3d vec);
	vec3d cross(vec3d pr);

	double norm();
	point to_point();

	double x, y, z;

	static bool collinear(vec3d a, vec3d b);

private:
	void init_coord(double s_x, double s_y, double s_z);
};

class edge {
 public:
	 edge();
	 edge(int sstart, int send);
	 bool operator == (const edge pr);
	 int number;

	 int start();
	 int end();


private:
	 pair<int, int> edge_node;
	 
};