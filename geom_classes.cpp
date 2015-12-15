#include "geom_classes.h"

point::point(){
	x = y = z = 0;
}

point::point(double sx, double sy, double sz){
	x = sx; y = sy; z = sz;
}

point::point(node nd){
	x = nd.x; y = nd.y; z = nd.z;
}

double& point::operator [] (int i){
	switch(i){
	case 0: return x;
	case 1: return y;
	case 2: return z;
	}
}

point& point::operator = (const point pr) {
	x = pr.x; y = pr.y; z=pr.z;
	return *this;
}

point& point::operator = (const node pr) {
	x = pr.x; y = pr.y; z=pr.z;
	return *this;
}

node::node(){
	x = y = z = 0;
	number = 0;
}

node::node(double sx, double sy, double sz, int num){
	x = sx; y = sy; z = sz;
	number = num;
}

double& node::operator [] (int i){
	switch(i){
	case 0: return x;
	case 1: return y;
	case 2: return z;
	}
}

node& node::operator = (const node pr) {
	x = pr.x; y = pr.y; z=pr.z; number=pr.number;
	return *this;
}

bool node::operator < (const node pr) {
	return number < pr.number ? true : false;
}

bool node::operator > (const node pr) {
	return  number > pr.number ? true : false;
}

ostream& operator << (ostream& os, node& out_n) {
	os << out_n.number << "\t" << out_n.x << "\t" << out_n.y << "\t" << out_n.z;
	return os;
}

istream& operator >> (istream& is, node& inp_n) {
	is >> inp_n.number >> inp_n.x >> inp_n.y >> inp_n.z;
	return is;
}

edge::edge() {
	edge_node.first = edge_node.second = 0;
}

edge::edge(int sstart, int send) {
	edge_node.first = sstart; edge_node.second = send;
}

bool edge::operator == (const edge ed) {
	if(edge_node == ed.edge_node)
		return true;
	else
		return false;
}

int edge::start() {
	return edge_node.first;
}

int edge::end() {
	return edge_node.second;
}

vec3d::vec3d() {
	x = y = z = 0;
}

vec3d::vec3d(point start, point end){
	x = end.x - start.x;
	y = end.y - start.y;
	z = end.z - start.z;
}

vec3d::vec3d(node start, node end){
	x = end.x - start.x;
	y = end.y - start.y;
	z = end.z - start.z;
}

vec3d::vec3d(double sx, double sy, double sz) {
	init_coord(sx, sy, sz);
}

vec3d::vec3d(node nd) {
	init_coord(nd.x, nd.y, nd.z);
}

void vec3d::init_coord(double s_x, double s_y, double s_z) {
	x = s_x; y = s_y; z = s_z;	
}

// Проверим, что вектора коллинеарны, т.е. линейно независимы
bool vec3d::collinear(vec3d a, vec3d b) {
	int cur = -1;
	bool found_cur = false;
	for(int i = 0; i < 3 && !found_cur; i++) {
		if (fabs(a[i]) > 1e-10) {
			cur = i;
			found_cur = true;
		}
	}
	if (cur == -1)
		return true;
	
	vec3d tmp_vec = (b[cur] / a[cur]) * a - b;

	bool collin[3];
	bool res = false;

	for(int i = 0; i < 3; i++) {
		collin[i] = fabs(tmp_vec[i]) > 1e-10;
		res = res || collin[i];
	}

	return !res;
}

double& vec3d::operator [] (const int i){
	switch(i){
	case 0: return x;
	case 1: return y;
	case 2: return z;
	}
}

double vec3d::operator * (vec3d pr) {
	return x*pr.x +y*pr.y + z*pr.z;
}

double vec3d::norm() {
	return sqrt(x*x + y*y + z*z);
}

point vec3d::to_point() {
	return point(x, y, z);
}

vec3d operator * (const double a, vec3d vec){
	return vec3d(a*vec.x, a*vec.y, a*vec.z);
}


vec3d operator * (const matrix(3) M, vec3d vec) {
	vec3d res(0, 0, 0);

	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			res[i] += M[i][j] * (vec[j]);
		}
	}

	return res;
}

vec3d vec3d::operator + (vec3d pr) {
	return vec3d(x+pr.x, y+pr.y, z+pr.z);
}

vec3d vec3d::operator - (vec3d pr) {
	return vec3d(x-pr.x, y-pr.y, z-pr.z);
}

bool vec3d::operator == (const vec3d pr) {
	return (x == pr.x) && (y == pr.y) && (z == pr.z);
}

bool vec3d::operator != (const vec3d pr) {
	return !this->operator==(pr);
}


vec3d vec3d::operator / (const double a) {
	return vec3d(x/a, y/a, z/a);
}

vec3d vec3d::cross(vec3d pr) {
	return vec3d(y*pr.z - z*pr.y, z*pr.x - x*pr.z, x*pr.y - y*pr.x);
}