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

	vector<s_vbfunc> f;

	s.get_basis(1, 1, f);


	return 0;
}