#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "geom_classes.h"
#include "L_cords_gen.h"
#include <functional>
#include <cstdio>
#include <vector>

typedef function<double(double, double, double)> func3d;
typedef function<vec3d(double, double, double)> vfunc3d;

const int gauss_points_tr = 3;
const int gauss_points_tet = 4;

class trelement;
class tetelement;

typedef double (trelement::*tr_bfunc)(double x, double y, double z);
typedef vec3d (trelement::*tr_vbfunc)(double x, double y, double z);
typedef double (tetelement::*tet_bfunc)(double x, double y, double z);
typedef vec3d (tetelement::*tet_vbfunc)(double x, double y, double z);


class sector {
public:

	 sector();
	 sector(vector<node> nodes_s, vector<dof_type> s_dofs);

	int& operator [] (int i); //�������� i-�� ������� �������

	static const int element_nodes = 2;

	double L2_diff(func3d f, vector<double>& q_loc);
private:
	vector<node> nodes;
	vector<int> dofs;

};

class trelement {
 public:
	 trelement();
	 trelement(vector<node> nodes_s, vector<dof_type> s_dofs);

	 int& operator [] (int i); //�������� i-�� ������� �������

	 double integral(func3d func); //���������� ��������� �� ������������ (�����, 3 �����) � ��������� �����������

	 int get_ph_area();
	 void set_ph_area(int sph_area);

	 dyn_matrix get_local_matrix(double mu);
	 vector<double> get_local_right_part(func3d rp_func);

	 double scalar_basis_v(int i, double x, double y, double z);

	 array<int, 3> loc_edge;

	 double integrate(func3d integ_func);//���������� ��������� �� ������������ � ���������� �����������

	 bool in_element(double x, double y, double z);
	 bool valid_for_tree_node(double x0, double x1, double y0, double y1, double z0, double z1);

	 vec3d basis_v(int i, double x, double y, double z);
	 node local_node(int i);
	 vec3d get_tau(int i);

	 vector<dof_type> get_dofs();

	 static const int element_nodes = 3;

	 double L2_diff(func3d f, vector<double>& q_loc);

	 double vector_jump_L2(vfunc3d f1, vfunc3d f2);
	 
 private:

	 void generate_L_cords();

	 unsigned int dofs_number;

	 point to_local_cord(point p_glob);
	 point to_global_cord(point p_loc);

	 vec3d to_global_cord(vec3d v_loc);

	 double lambda(int l_i, point p_loc);

	 vec3d normal_vector;

	 vector<tr_bfunc> scalar_basis;
	 vector<tr_vbfunc> scalar_basis_grad;

	 vector<dof_type> dofs;

	 matrix(3) D_matrix, L_cord_matrix; //L - ������� L-���������, D - � ���������
	 double det_D; //������������� ������� L-���������

	 void init_cords(); //���������� ��������� ���������

	  array<node, 3> node_array;
	  array<int, 3> edge_array;


	 /* ��������� �����
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

	  matrix(3) transition_matrix; //������� �������� � ��������� ����������

	  point transform(point pr); //��������� ���������� ���������� � ���������

	  int ph_area; //���������� �������
	  point trpoint[3]; //��������� ���������� ����� ������������
	  point gauss_points[gauss_points_tr]; //����� ��� �������������� �� ������ (� ��������� ������� ���������)
	  point gauss_points_global[gauss_points_tr]; //����� ��� �������������� �� ������ (� ��������� ������� ���������)
	  double gauss_weights[gauss_points_tr];
	  double jacobian; //������� ��� ��������� ���������

	  array<vec3d, 4> tau;

	  double basis_1(double x, double y, double z);
	  double basis_2(double x, double y, double z);
	  double basis_3(double x, double y, double z);

	  vec3d grad_basis_1(double x, double y, double z);
	  vec3d grad_basis_2(double x, double y, double z);
	  vec3d grad_basis_3(double x, double y, double z);
};


class tetelement {
 public:
	 tetelement();
	 tetelement(vector<node> nodes_s, vector<dof_type> s_dofs);

	 int& operator [] (int i); //�������� i-�� ��������� �����

	 dyn_matrix get_local_matrix(double mu);
	 vector<double> get_local_right_part(func3d rp_func);

	 double integrate(func3d integ_func);

	 int get_ph_area();
	 void set_ph_area(int sph_area);

	 node get_local_node(int i);
	 point get_center();
	 double scalar_basis_v(int i, double x, double y, double z);
	 vec3d scalar_basis_grad_v(int i, double x, double y, double z);
	 double scalar_basis_div_grad_v(int i, double x, double y, double z);

	 bool in_element(double x, double y, double z);
	 bool valid_for_tree_node(double x0, double x1, double y0, double y1, double z0, double z1);

	 vector<dof_type> get_dofs();

	 static const int element_nodes = 4;

	 double L2_diff(func3d f, vector<double>& q_loc, vector<double>& q_virtual);

 private:

	 vector<tet_bfunc> scalar_basis;
	 vector<tet_vbfunc> scalar_basis_grad;
	 vector<tet_bfunc> scalar_basis_div_grad;
	 vector<dof_type> dofs;
	 unsigned int dofs_number;

	 void init_coords();
	 void generate_L_cords();
	 void calculate_M_matrix();

	 double lambda(int i, point p_glob);	//�������� i-� L-���������� � �����
	 vec3d grad_lambda(int i);				//�������� i-� L-����������

	 array<node,4> node_array;
	 array<int,6> edge_array;

	 int ph_area; //���������� �������


	 matrix(4) D_matrix, L_cord_matrix; //L - ������� L-���������, D - � ���������
	 //i - ������ �� i-� ����������, ���� ������������ x,y,z,1
	 double det_D; //������������� ������� L-���������
	 matrix(12) M_matrix;

	  point gauss_points[gauss_points_tet]; //����� ��� �������������� �� ������
	  double gauss_weights[gauss_points_tet];
	  double jacobian; //������� ��� ��������� ���������

	  //��� ���������� ������
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