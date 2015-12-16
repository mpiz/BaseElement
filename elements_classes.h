#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "add_geom_classes.h"
#include "L_cords_gen.h"
#include <functional>
#include <cstdio>
#include <vector>

typedef function<double(double, double, double)> func3d;
typedef function<vec3d(double, double, double)> vfunc3d;

const int gauss_points_sec = 3;
const int gauss_points_tr = 3;
const int gauss_points_tet = 4;

class simple_element;
class sector;
class trelement;
class tetelement;


typedef double (trelement::*tr_bfunc)(double x, double y, double z);
typedef vec3d (trelement::*tr_vbfunc)(double x, double y, double z);
typedef double (tetelement::*tet_bfunc)(double x, double y, double z);
typedef vec3d (tetelement::*tet_vbfunc)(double x, double y, double z);

class simple_element {

public:
	virtual vfunc3d get_vector_basis_dof(size_t dof_i);
	virtual vfunc3d get_vector_basis(dof_type order, dof_type num = 0);


	virtual vfunc3d get_vector_right_part_dof(size_t dof_i);
	virtual vfunc3d get_vector_right_part(dof_type order, dof_type num = 0);

	virtual dyn_matrix get_local_matrix(double mu);
	virtual vector<double> get_local_right_part(func3d rp_func);
	virtual vector<double> get_local_right_part(vfunc3d rp_func);
	virtual double integrate(func3d func);	// ¬ычисление интеграла по элементу


	void add_dof(dof_info d);
	void prepare_gauss(int gn);
	vector<dof_info> get_dofs();
	vector<dof_type> get_dofs_num();


protected:
	
	vector<dof_info> dofs;
	unsigned int dofs_number;


	size_t gauss_points_n;
	double jacobian;
	vector<point> gauss_points_global; //точки дл€ интегрировани€ по √ауссу (в локальной системе координат)
	vector<double> gauss_weights;
};


class sector : public simple_element {
public:


	sector();
	sector(const vector<node>& nodes_s, const plane& plane_s);
	sector(vector<node> nodes_s, vector<dof_info> s_dofs);

	int& operator [] (int i); //получить i-ую степень свободы

	static const int element_nodes = 2;

	double L2_diff(func3d f, vector<double>& q_loc);

	vfunc3d get_vector_basis_dof(size_t dof_i);
	vfunc3d get_vector_basis(dof_type order, dof_type num);
	func3d get_vector_basis_dof_tau(size_t dof_i);		// “ангенсальна€ компонента функции get_vector_basis_dof(dof_i)

	vfunc3d get_vector_right_part_dof(size_t dof_i);
	vfunc3d get_vector_right_part(dof_type order, dof_type num = 0);


	static dof_type get_dof_n(dof_type order, dof_type num);

	dyn_matrix get_local_matrix(double mu);
	vector<double> get_local_right_part(vfunc3d rp_func);

	bool in_element(double x, double y, double z);

	void for_point_on_element(function<void(double, double, double)> func);

private:
	vector<node> nodes;
	plane sector_plane;	// ѕлоскость, в которой лежит отрезок

	vec3d direction;
	double length;

	vec3d normal_in_plane;	// Ќормальный вектор в плоскости sector_plane

	double get_t(double x, double y, double z);
	point get_point(double t);

	vec3d vbasis_1_1(double x, double y, double z);
	vec3d vbasis_1_2(double x, double y, double z);

	double k_sq;

	void init_coords();




};

class trelement : public simple_element {
 public:
	 trelement();
	 trelement(vector<node> nodes_s, vector<dof_info> s_dofs);

	 int& operator [] (int i); //получить i-ую степень свободы

	 double integral(func3d func); //вычисление интеграла по треугольнику (√аусс, 3 точки) в локальныъ координатах

	 int get_ph_area();
	 void set_ph_area(int sph_area);

	 dyn_matrix get_local_matrix(double mu);
	 vector<double> get_local_right_part(func3d rp_func);
	 vector<double> get_local_right_part(vfunc3d rp_func);

	 double scalar_basis_v(int i, double x, double y, double z);

	 array<int, 3> loc_edge;

	 bool in_element(double x, double y, double z);
	 bool valid_for_tree_node(double x0, double x1, double y0, double y1, double z0, double z1);

	 vec3d basis_v(int i, double x, double y, double z);
	 node local_node(int i);
	 vec3d get_tau(int i);

	 static const int element_nodes = 3;

	 double L2_diff(func3d f, vector<double>& q_loc);

	 double vector_jump_L2(vfunc3d f1, vfunc3d f2);
	 vfunc3d get_vector_basis_dof(size_t dof_i);

	 static dof_type get_dof_n(dof_type order, dof_type num);
	 
 private:

	 plane tr_plane;

	 void generate_L_cords();

	 point to_local_cord(point p_glob);
	 point to_global_cord(point p_loc);

	 vec3d to_global_cord(vec3d v_loc);

	 double lambda(int l_i, point p_loc);

	 vec3d normal_vector;

	 vector<tr_bfunc> scalar_basis;
	 vector<tr_vbfunc> scalar_basis_grad;

	 matrix(3) D_matrix, L_cord_matrix; //L - мтарица L-координат, D - еЄ обратнна€
	 double det_D; //опеределитель матрицы L-координат

	 void init_cords(); //построение локальных координат

	  array<node, 3> node_array;
	  array<int, 3> edge_array;


	 /* локальный базис
	   w1 = l1 * grad(l2) - l2 * grad(l1);
	  w2 = l1 * grad(l3) - l3 * grad(l2);
	  w3 = l2 * grad(l3) - l3 * grad(l2);

	  w4 = l1 * grad(l2) + l2 * grad(l1);
	  w5 = l1 * grad(l3) + l3 * grad(l2);
	  w6 = l2 * grad(l3) + l3 * grad(l2);
	  */

	  vec3d gradphi1(point p_loc);
	  vec3d gradphi2(point p_loc);
	  vec3d gradphi3(point p_loc);

	  vec3d grad_lambda(int i);

	  matrix(3) transition_matrix; //матрица перехода в локальные координаты


	  int ph_area; //физическа€ область
	  array<point,3> trpoint; //локальные координаты точек треугольника
	  vector<point> gauss_points; //точки дл€ интегрировани€ по √ауссу (в локальной системе координат)
	  double jacobian; //€кобиан дл€ вычилени€ интеграла


	  double basis_1(double x, double y, double z);
	  double basis_2(double x, double y, double z);
	  double basis_3(double x, double y, double z);

	  vec3d grad_basis_1(double x, double y, double z);
	  vec3d grad_basis_2(double x, double y, double z);
	  vec3d grad_basis_3(double x, double y, double z);

	  vector<vfunc3d> vector_basis;

	  vfunc3d get_vector_basis_for_dof(dof_type order, dof_type num, dof_type n1, dof_type n2);
};


class tetelement : public simple_element {
 public:
	 tetelement();
	 tetelement(vector<node> nodes_s, vector<dof_info> s_dofs);

	 int& operator [] (int i); //получить i-ое локальное ребро

	 dyn_matrix get_local_matrix(double mu);
	 vector<double> get_local_right_part(func3d rp_func);

	 int get_ph_area();
	 void set_ph_area(int sph_area);

	 node get_local_node(int i);
	 point get_center();
	 double scalar_basis_v(int i, double x, double y, double z);
	 vec3d scalar_basis_grad_v(int i, double x, double y, double z);
	 double scalar_basis_div_grad_v(int i, double x, double y, double z);

	 bool in_element(double x, double y, double z);
	 bool valid_for_tree_node(double x0, double x1, double y0, double y1, double z0, double z1);

	 static const int element_nodes = 4;

	 double L2_diff(func3d f, vector<double>& q_loc, vector<double>& q_virtual);

 private:

	 vector<tet_bfunc> scalar_basis;
	 vector<tet_vbfunc> scalar_basis_grad;
	 vector<tet_bfunc> scalar_basis_div_grad;

	 void init_coords();
	 void generate_L_cords();
	 void calculate_M_matrix();

	 double lambda(int i, point p_glob);	//значение i-й L-координаты в точке
	 vec3d grad_lambda(int i);				//градиент i-й L-координаты

	 array<node,4> node_array;
	 array<int,6> edge_array;

	 int ph_area; //физическа€ область


	 matrix(4) D_matrix, L_cord_matrix; //L - матрица L-координат, D - еЄ обратнна€
	 //i - строка за i-ю координану, коэф соответсенно x,y,z,1
	 double det_D; //опеределитель матрицы L-координат
	 matrix(12) M_matrix;

	  //дл€ построени€ дерева
	  array<double, 3> ch_points[5];
	  double edges_a[6][3], edges_b[6][3];

	  double basis_1(double x, double y, double z);
	  double basis_2(double x, double y, double z);
	  double basis_3(double x, double y, double z);
	  double basis_4(double x, double y, double z);


	  vec3d grad_basis_1(double x, double y, double z);
	  vec3d grad_basis_2(double x, double y, double z);
	  vec3d grad_basis_3(double x, double y, double z);
	  vec3d grad_basis_4(double x, double y, double z);

	  double div_grad_basis_1(double x, double y, double z);
	  double div_grad_basis_2(double x, double y, double z);
	  double div_grad_basis_3(double x, double y, double z);
	  double div_grad_basis_4(double x, double y, double z);
};