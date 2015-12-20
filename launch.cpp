#include "VirtualFace_Hrot.h"

int main() {

	vector<node> nodes;
	nodes.push_back(node(0, 0, 0, 0));
	nodes.push_back(node(1, 0, 0, 1));
	nodes.push_back(node(0, 1, 0, 2));

/*	nodes.push_back(node(0.2, 1, 0, 3));
	nodes.push_back(node(0.8, 1, 0, 4));

	nodes.push_back(node(1, 0.8, 0, 5));
	nodes.push_back(node(1, 0.2, 0, 6));

	nodes.push_back(node(0.8, 0, 0, 7));*/

	VirtualFace_Hrot face(nodes, 1, 0);

	face.input_mesh("FaceTriangle.dat");
	face.input_bound("FaceTriangleBound.dat");

	face.calculate();

	face.test_print_local_basis("basis");
	face.test_calc_points(2);
	//face.test_func_info(0);
	//face.test_func_info(1);
	face.test_func_info(2);

	return 0;
}