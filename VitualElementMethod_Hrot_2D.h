#pragma once
#include "VirtualFace_Hrot.h"

const int virtual_method_order = 1;

class VirtualElementMethod_Hrot_2D : public BaseElement<VirtualFace_Hrot> {
public:
	VirtualElementMethod_Hrot_2D();

	void calculate();
	void input_mesh(string file_name);
	void input_bound(string file_name);



	vector<dof_info> calc_element_dofs(vector<node>& el_nodes);

	vec3d vector_basis_val(dof_type basis_i, double x, double y, double z);
	VirtualFace_Hrot* find_element(point pn);

private:

	int local_edges_n;
	vector<sector*> local_edges;

	map_node_to_edge nodes_to_edges;

};