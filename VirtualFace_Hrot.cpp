#include "VirtualFace_Hrot.h"


VirtualFace_Hrot::VirtualFace_Hrot() {
}

VirtualFace_Hrot::VirtualFace_Hrot(const vector<node>& nodes_s, const vector<dof_type>& dofs_s, size_t order, dof_type num = 0) {

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

	//—формируем рЄбра
	for (size_t node_i = 0; node_i < nodes_n; node_i) {
		size_t node_next = node_i%nodes_n;
		sector_nodes[0] = nodes[node_i];
		sector_nodes[1] = nodes[node_next];

		edges.push_back(sector(sector_nodes, face_plane));
	}

}

VirtualFace_Hrot::~VirtualFace_Hrot() {

}
