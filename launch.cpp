#include "VitualElementMethod_Hrot_2D.h"

int main() {

	vector<node> nodes;
	//nodes.push_back(node(0, 0, 0, 0));
	//nodes.push_back(node(1, 0, 0, 1));
	//nodes.push_back(node(0, 1, 0, 2));

	/*nodes.push_back(node(0.2, 0, 0, 0));

	nodes.push_back(node(0, 0.2, 0, 1));
	nodes.push_back(node(0, 0.8, 0, 2));

	nodes.push_back(node(0.2, 1, 0, 3));

	nodes.push_back(node(0.6, 0.6, 0, 4));

	nodes.push_back(node(0.8, 1, 0, 5));

	nodes.push_back(node(1, 0.8, 0, 6));
	nodes.push_back(node(1, 0.2, 0, 7));

	nodes.push_back(node(0.8, 0, 0, 8));

	map_node_to_edge gg;
	gg[pair<size_t, size_t>(0, 1)] = 0;
	gg[pair<size_t, size_t>(1, 2)] = 0;
	gg[pair<size_t, size_t>(2, 3)] = 0;
	gg[pair<size_t, size_t>(3, 4)] = 0;
	gg[pair<size_t, size_t>(4, 5)] = 0;
	gg[pair<size_t, size_t>(5, 6)] = 0;
	gg[pair<size_t, size_t>(6, 7)] = 0;
	gg[pair<size_t, size_t>(7, 8)] = 0;
	gg[pair<size_t, size_t>(0, 8)] = 0;

	VirtualFace_Hrot face(nodes, gg, 1, 0);

	face.input_mesh("Face2.dat");
	face.input_bound("FaceBound2.dat");

	face.calculate();
	face.test_print_local_basis("test1func.txt");
	face.test_calc_points(2);
	for(int i = 0; i < 9; i++)
		face.test_func_info(i);*/

	VirtualElementMethod_Hrot_2D method;

	method.input_mesh("virtual_mesh.txt");
	method.calculate();

	return 0;
}