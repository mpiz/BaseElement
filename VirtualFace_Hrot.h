#pragma once

#include "BaseElement.h"
#include "VirtualEdge_Hrot.h"



class VirtualFace_Hrot : public BaseElement<trelement> {
public:
	VirtualFace_Hrot();
	VirtualFace_Hrot(const vector<node>& nodes_s, dof_type order, dof_type num = 0);

	void calculate();
	void input_mesh(string file_name);
	void input_bound(string file_name);

	void test_calc_points(dof_type dof_i);


	~VirtualFace_Hrot();

private:

	plane face_plane;
	size_t edges_n;	// Количество ребёр

	vector<sector*> edges;

	vector<vfunc3d> bound_functions;	// Функции для учёта краевых условий
	vector<vfunc3d> right_part_functions; // Функции правой части

	vector<sector> bound_elements;


	VirtualEdge_Hrot bound_edge;

	double get_bound_value(dof_type basis_i, dof_type cur_dof);

	dof_type get_bound_dof(dof_type face_dof);

	vector<double*> bound_solutions;


};

