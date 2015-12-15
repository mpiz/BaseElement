#include "VirtualFace_Hrot.h"

int main() {

	vector<node> nodes;
	nodes.push_back(node(0.2, 0, 0, 0));
	nodes.push_back(node(0, 0.2, 0, 1));
	nodes.push_back(node(0, 0.8, 0, 2));

/*	nodes.push_back(node(0.2, 1, 0, 3));
	nodes.push_back(node(0.8, 1, 0, 4));

	nodes.push_back(node(1, 0.8, 0, 5));
	nodes.push_back(node(1, 0.2, 0, 6));

	nodes.push_back(node(0.8, 0, 0, 7));*/

	VirtualFace_Hrot face(nodes, 1, 1);

	face.input_mesh("Face.dat");
	face.input_bound("FaceBound.dat");

	face.calculate();

	face.test_calc_points(0);

	return 0;
}