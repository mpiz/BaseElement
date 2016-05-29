#include "VitualElementMethod_Hrot_2D.h"

int main() {

	vector<node> trnodes;
	vector<dof_info> trdofs;
	
	trnodes.push_back(node(0, 0, 0, 0));
	trnodes.push_back(node(0.2, 0, 0, 1));
	trnodes.push_back(node(0, 0.2, 0, 2));

	dof_info di[3];
	di[0].number = 0;
	di[0].order = 1;
	di[0].num = 0;
	di[0].add_geom(0);
	di[0].add_geom(1);

	di[1].number = 1;
	di[1].order = 1;
	di[1].num = 0;
	di[1].add_geom(0);
	di[1].add_geom(2);

	di[2].number = 2;
	di[2].order = 1;
	di[2].num = 0;
	di[2].add_geom(1);
	di[2].add_geom(2);

	trdofs.push_back(di[0]);
	trdofs.push_back(di[1]);
	trdofs.push_back(di[2]);

	trelement tr(trnodes, trdofs);

	auto A_loc = tr.get_local_matrix(0);

	VirtualElementMethod_Hrot_2D method;

	cout << "Input mesh\n";
	method.input_mesh("virtual_mesh.txt");
	cout << "Input bound\n";
	method.input_bound("virtual_mesh_bound.txt");
	method.calculate();
	cout << "Printing solution\n";
	method.test_print_solution("virtual_solution.txt");

	cout << "Calculating error\n";
	double err = method.calc_error();
	FILE* outp = fopen("error.txt", "w");
	fprintf(outp, "%.5e", err);
	fclose(outp);

	return 0;
}