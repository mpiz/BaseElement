#include "VirtualFace_Hrot.h"


VirtualFace_Hrot::VirtualFace_Hrot() : BaseElement() {

}

void VirtualFace_Hrot::input_mesh(string file_name) {
	BaseElement::input_mesh(file_name, 203);
}

void VirtualFace_Hrot::input_bound(string file_name) {
	bound_edge.input_mesh(file_name);
}

vector<dof_info> VirtualFace_Hrot::calc_element_dofs(vector<node>& el_nodes) {
	return calc_element_dofs_edge(el_nodes, sector::get_dof_n(method_order, method_num));
}

VirtualFace_Hrot::VirtualFace_Hrot(const vector<node>& nodes_s, size_t order, dof_type num) {

	nodes = nodes_s;
	nodes_n = nodes.size();
	method_order = order;
	method_num = num;

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

			dof_info add_dof(counters::get_next_dof(), order, num, dof_i);
			add_dof.add_geom(nodes[node_i].number);
			add_dof.add_geom(nodes[node_next].number);

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
	bound_edge.get_solutions(bound_solutions);

	auto firts_bound_geom = bound_edge.get_bound_funcs();
	for(auto& fb_g : firts_bound_geom) {
		first_bound.push_back(edge_type_dofs[fb_g].number);
	}

	sort(first_bound.begin(), first_bound.end());

	generate_port();
	generate_matrix_with_out_bound(right_part_functions);

	generate_matrix_first_bound([&](dof_type basis_i, dof_type cur_dof)->double {
		return get_bound_value(basis_i, cur_dof);
	});

#ifdef DEBUGOUTP
	print_full_matrix("matrix.txt");

	for(size_t i = 0; i < dofs_n; i++) {
		string s;
		stringstream ss;
		ss << "rp_" << i << ".txt";
		ss >> s;
		print_right_part(i, s);
	}
#endif
	
	solve_SLAE();

}

dof_type VirtualFace_Hrot::get_bound_dof(dof_type face_dof) {
	auto tup = edge_type_dofs_revers[face_dof];
	auto res = bound_edge.get_dof_num(tup);
	return res;
}

double VirtualFace_Hrot::get_bound_value(dof_type basis_i, dof_type cur_dof) {
	auto bound_dof = get_bound_dof(cur_dof);
	return bound_solutions[basis_i][bound_dof];
}


vec3d VirtualFace_Hrot::vector_basis_val(dof_type basis_i, double x, double y, double z){
	point pn(x, y, z);
	auto el = find_element(pn);
	if (el == nullptr)
		return vec3d(0, 0, 0);

	auto el_dofs = el->get_dofs();
	auto el_dofs_n = el_dofs.size();

	vec3d res(0, 0, 0);
	double* basis_sol = solutions[basis_i];

	for(size_t k = 0; k < el_dofs_n; k++) {
		vec3d basis_val = el->get_vector_basis_dof(k)(x, y, z);
		res = res + basis_sol[el_dofs[k].number] * basis_val;
	}
	return res;
}

trelement* VirtualFace_Hrot::find_element(point pn) {
	for(int el_i = 0; el_i < elements_n; el_i++) {
		if (elements[el_i]->in_element(pn.x, pn.y, pn.z))
			return elements[el_i];
	}
	return nullptr;
}

void VirtualFace_Hrot::test_calc_points(dof_type dof_i) {
	bound_edge.print_full_matrix("test_matrix.txt");
	bound_edge.print_right_part(dof_i, "test_rp.txt");
	bound_edge.test_calc_points(dof_i);

}

void VirtualFace_Hrot::test_func_info(dof_type dof_i) {
	string file_name;
	stringstream ss;
	ss << "test_2_" << dof_i << ".txt";
	ss >> file_name;


	ofstream outpfile(file_name.c_str());

	double y = 0;
	double h = 0.02;

	double* q = solutions[dof_i];

	outpfile << "VARIABLES = \"x\" \"y\" \"v1\" \"v2\" \"t1\" \"t2\" \"t3\" \n";  
	
	while(y < 1) {
		double x = 0;
		while (x <= 1) {
			vec3d val = vector_basis_val(dof_i, x, y, 0);
			double tau1 = val * vec3d(1, 0, 0);
			double tau2 = val * vec3d(-1, 1, 0) / vec3d(-1, 1, 0).norm();
			double tau3 = val * vec3d(0, 1, 0);
			outpfile << x << " " << y << " " << val.x << " " << val.y << " " << tau1 << " " << tau2 << " " << tau3 << endl; 
			x += h;
		}
		y += h;
	}

	outpfile.close();

}

void VirtualFace_Hrot::test_print_local_basis(string file_name) {
	ofstream outp(file_name.c_str());

	double y = 0;
	double h = 0.02;
	outp << "VARIABLES = \"x\" \"y\" ";  
	for(auto i = 0; i < local_dof_n; i++) {
		outp << "\"vbasis_" << i << "_x\" ";
		outp << "\"vbasis_" << i << "_y\" ";
	}
	outp << endl;

	while(y < 1) {
		double x = 0;
		while (x <= 1) {
			outp << x << " " << y << " ";
			auto el_it = find_element(point(x, y, 0));
			for(auto& dof_it : local_dofs) {
				vec3d res1(0, 0, 0);
				auto el_dof = el_it->get_dofs_num();
					for(size_t dof_i = 0; dof_i < el_dof.size(); dof_i++) {
						if (el_dof[dof_i] == dof_it.number) {
							res1 = el_it->get_vector_basis_dof(dof_i)(x, y, 0);
							break;
						}
					}
				outp << res1.x << " " << res1.y << " ";
			}
			outp << endl;
			x += h;
		}
		y += h;
	}

	outp.close();

}

VirtualFace_Hrot::~VirtualFace_Hrot() {
	for(auto& edge : edges)
		delete edge;
}