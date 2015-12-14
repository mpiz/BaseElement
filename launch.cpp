#include "elements_classes.h"

int main() {

	array<node, 3> nodes;
	nodes[0] = node(0, 0, 0);
	nodes[1] = node(1, 0, 0);
	nodes[2] = node(0, 1, 0);

	plane pl(nodes);

	vector<node> s_nodes;
	s_nodes.push_back(node(0, 0, 0));
	s_nodes.push_back(node(1, 1, 0));

	sector s(s_nodes, pl);

	vector<vbfunc> f;
	f.reserve(2);

	f.push_back(s.get_vector_basis(1, 1));
	f.push_back(s.get_vector_basis(1, 2));

	auto res = f[0](0.5, 0.5, 0);
	auto res2 = f[1](0.5, 0.5, 0);


	auto val = res * vec3d(1, 1, 0) / vec3d(1, 1, 0).norm();
	auto val2 = res2 * vec3d(1, 1, 0) / vec3d(1, 1, 0).norm();

	return 0;
}