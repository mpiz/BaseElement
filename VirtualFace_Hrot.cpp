#include "VirtualFace_Hrot.h"


VirtualFace_Hrot::VirtualFace_Hrot() : BaseElement() {

}

void VirtualFace_Hrot::input_mesh(string file_name) {
	BaseElement::input_mesh(file_name, 203);
}

void VirtualFace_Hrot::input_bound(string file_name) {
	bound_edge.input_mesh(file_name);
}

VirtualFace_Hrot::VirtualFace_Hrot(const vector<node>& nodes_s, size_t order, dof_type num) {

	nodes = nodes_s;
	nodes_n = nodes.size();
	element_order = order;

	if (nodes_n < 3) {
		throw "VirtualFace_Hrot::VirtualFace_Hrot - to construct a face atleast 3 points needed";
	}

	// —формируем грань, как плоскость
	array<node, 3> plane_nodes;
	bool found_plane = false;
	for (size_t i = 0; i < nodes_n && !found_plane; i++) {
		plane_nodes[0] = nodes[i];
		for (size_t j = i + 1; j < nodes_n && !found_plane; j++) {
			plane_nodes[1] = nodes[j];
			for (size_t k = j + 1; k < nodes_n && !found_plane; k++) {
				plane_nodes[2] = nodes[k];
				found_plane = !plane::is_on_line(plane_nodes);
			}
				
		}
	}

	if (!found_plane) {
		throw "VirtualFace_Hrot::VirtualFace_Hrot - there is no plane, there is a line";
	}

	face_plane = plane(plane_nodes);

	edges_n = nodes_n;
	edges.reserve(edges_n);
	vector<node> sector_nodes;
	sector_nodes.resize(2);

	dofs_n = 0;

	//—формируем рЄбра
	for (size_t node_i = 0; node_i < nodes_n; node_i++) {
		size_t node_next = (node_i + 1)%nodes_n;
		sector_nodes[0] = nodes[node_i];
		sector_nodes[1] = nodes[node_next];

		sector* add_sector = new sector(sector_nodes, face_plane);

		auto sec_dof_n = sector::get_dof_n(order, num);
		for(size_t dof_i = 0; dof_i < sec_dof_n; dof_i++) {
			dof_type add_dof = counters::get_next_dof();
			dofs.push_back(add_dof);
			add_sector->add_dof(add_dof);

			bound_functions.push_back(add_sector->get_vector_basis_dof(dof_i));
			right_part_functions.push_back(add_sector->get_vector_right_part_dof(dof_i));
		}

		edges.push_back(add_sector);

	}

	dofs_n = dofs.size();
	bound_edge = VirtualEdge_Hrot(order, num);

}

void VirtualFace_Hrot::calculate() {

	bound_edge.set_right_parts(bound_functions);
	bound_edge.calculate();

}

void VirtualFace_Hrot::test_calc_points(dof_type dof_i) {
	bound_edge.print_full_matrix("test_matrix.txt");
	bound_edge.print_right_part(dof_i, "test_rp.txt");
	bound_edge.test_calc_points(dof_i);


}

VirtualFace_Hrot::~VirtualFace_Hrot() {
	for(auto& edge : edges)
		delete edge;
}
