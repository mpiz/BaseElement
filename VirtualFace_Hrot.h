#pragma once

#include "BaseElement.h"
#include "VirtualEdge_Hrot.h"

typedef map<pair<size_t, size_t>, size_t> map_node_to_edge;


class VirtualFace_Hrot : public BaseElement<trelement> {
public:
	VirtualFace_Hrot();
	VirtualFace_Hrot(const vector<node>& nodes_s, map_node_to_edge& node_to_edge, dof_type order, dof_type num = 0);

	void calculate();
	void input_mesh(string file_name);
	void input_bound(string file_name);

	void input_mesh_from_params();

	vector<dof_info> calc_element_dofs(vector<node>& el_nodes);

	vec3d vector_basis_val(dof_type basis_i, double x, double y, double z);
	trelement* find_element(point pn);

	vector<dof_type> get_dofs_num();
	
	void test_calc_points(dof_type dof_i);
	void test_func_info(dof_type dof_i);
	void test_print_local_basis(string file_name);

	void set_mesh_files(string main_mesh, string bound_mesh);

	dyn_matrix get_local_matrix(double k_sq);

	vector<dof_info> get_dofs();

	vector<double> get_local_right_part(vfunc3d rp_func);
	vec3d get_vector_basis_dof(dof_type& dof_i, double x, double y, double z);

	double integrate(func3d func);

	~VirtualFace_Hrot();

private:

	plane face_plane;
	size_t edges_n;	// Количество ребёр

	vector<sector*> edges;

	vector<vfunc3d> bound_functions;	// Функции для учёта краевых условий
	vector<vfunc3d> right_part_functions; // Функции правой части

	vector<sector*> bound_elements;


	VirtualEdge_Hrot bound_edge;

	double get_bound_value(dof_type basis_i, dof_type cur_dof);

	dof_type get_bound_dof(dof_type face_dof);

	vector<double*> bound_solutions;

	string main_mesh_file_name, bound_mesh_file_name;


};

