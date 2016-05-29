#include "VirtualFace_Hrot.h"

const double mu0_inv = 1.0/(16 * atan(1.0) * 1e-7);


VirtualFace_Hrot::VirtualFace_Hrot() : BaseElement() {

}

void VirtualFace_Hrot::input_mesh(string file_name) {
	BaseElement::input_mesh(file_name, 203);
}

void VirtualFace_Hrot::input_bound(string file_name) {
	bound_edge.input_mesh(file_name);
}

void VirtualFace_Hrot::input_mesh_from_params() {
	input_mesh(main_mesh_file_name);
	input_bound(bound_mesh_file_name);
}

vector<dof_info> VirtualFace_Hrot::calc_element_dofs(vector<node>& el_nodes) {
	return calc_element_dofs_edge(el_nodes, sector::get_dof_n(method_order, method_num));
}

VirtualFace_Hrot::VirtualFace_Hrot(const vector<node>& nodes_s, map_node_to_edge& node_to_edge, dof_type order, dof_type num) {

	nodes = nodes_s;
	nodes_n = nodes.size();
	method_order = order;
	method_num = num;

	if (nodes_n < 3) {
		throw "VirtualFace_Hrot::VirtualFace_Hrot - to construct a face atleast 3 points needed";
	}

	// Сформируем грань, как плоскость
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

	//Сформируем рёбра
	for (size_t node_i = 0; node_i < nodes_n; node_i++) {
		size_t node_next = (node_i + 1)%nodes_n;
		sector_nodes[0] = nodes[node_i];
		sector_nodes[1] = nodes[node_next];

		sector* add_sector = new sector(sector_nodes, face_plane);

		auto global_n1 = nodes[node_i].number;
		auto global_n2 = nodes[node_next].number;
		if (global_n1 > global_n2)
			swap(global_n1, global_n2);

		auto global_dof = node_to_edge[pair<size_t, size_t>(global_n1, global_n2)];

		// Вводим упрощение, что sec_dof_n = 1
		auto sec_dof_n = sector::get_dof_n(order, num);
		if (sec_dof_n != 1)
			throw "VirtualFace_Hrot::VirtualFace_Hrot - only 1 node for each edge can be used";

		for(size_t dof_i = 0; dof_i < sec_dof_n; dof_i++) {

			dof_info add_dof(counters::get_next_dof(this), order, num, dof_i);
			add_dof.add_geom(nodes[node_i].number);
			add_dof.add_geom(nodes[node_next].number);

			dofs.push_back(add_dof);
			add_sector->add_dof(add_dof);

			loc_to_glob[add_dof.number] = global_dof;
			glob_to_loc[global_dof] = add_dof.number;

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


dyn_matrix VirtualFace_Hrot::get_local_matrix(double k_sq) {
	dyn_matrix A_loc;
	A_loc.resize(dofs_n);

	for (int i = 0; i < dofs_n; i++) {
		A_loc[i].resize(dofs_n);
		for (int j = 0; j < dofs_n; j++)
			A_loc[i][j] = 0;
	}
	
	//Пробежимcя по элементам, посчитаем внутри каждого квадратичную форму и сложим
	for (auto el_i : elements) {
		for (int i = 0; i < dofs_n; i++) {
			for (int j = 0; j <= i; j++) {

				auto el_dofs = el_i->get_dofs_num();
				int el_dofs_n = el_dofs.size();

				double el_res = el_i->integrate([&](double x, double y, double z)->double {
					vec3d rot_i(0, 0, 0), rot_j(0, 0, 0), val_i(0, 0, 0), val_j(0, 0, 0);
					for (int tr_dof_i = 0; tr_dof_i < el_dofs_n; tr_dof_i++) {
						size_t el_i_dof = el_dofs[tr_dof_i];
						double q_i = solutions[i][el_i_dof];
						double q_j = solutions[j][el_i_dof];
						vec3d el_basis = el_i->get_vector_basis_dof(tr_dof_i)(x, y, z);
						vec3d el_rot = el_i->get_vector_basis_rot_dof(tr_dof_i)(x, y, z);

						rot_i = rot_i + q_i * el_rot;
						rot_j = rot_j + q_j * el_rot;

						val_i = val_i + q_i * el_basis;
						val_j = val_j + q_j * el_basis;	
					}
					double add_el_value = mu0_inv * rot_i * rot_j + k_sq * val_i * val_j;
					return add_el_value;
				});

				A_loc[i][j] += el_res;
				if(i != j)
					A_loc[j][i] += el_res;
			}
		}

	}

	return A_loc;

}

vector<double> VirtualFace_Hrot::get_local_right_part(vfunc3d rp_func) {
	vector<double> b_loc;
	b_loc.resize(dofs_n);
	for (int i = 0; i < dofs_n; i++)
		b_loc[i] = 0;

	//Пробежимcя по элементам, посчитаем внутри каждого линейную форму и сложим
	for (auto el_i : elements) {
		for (int i = 0; i < dofs_n; i++) {

			auto el_dofs = el_i->get_dofs_num();
			int el_dofs_n = el_dofs.size();

			double el_res = el_i->integrate([&](double x, double y, double z)->double {
				vec3d val_i(0, 0, 0);
				for (int tr_dof_i = 0; tr_dof_i < el_dofs_n; tr_dof_i++) {
					size_t el_i_dof = el_dofs[tr_dof_i];
					double q_i = solutions[i][el_i_dof];
					vec3d el_basis = el_i->get_vector_basis_dof(tr_dof_i)(x, y, z);

					val_i = val_i + q_i * el_basis;

				}
				vec3d val_f = rp_func(x, y, z);
				double add_el_value = val_i * val_f;
				return add_el_value;
			});

			b_loc[i] += el_res;
		}

	}

	return b_loc;

}

vector<dof_info> VirtualFace_Hrot::get_dofs() {
	return dofs;
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
	for(auto i = 0; i < dofs_n; i++) {
		outp << "\"vbasis_" << i << "_x\" ";
		outp << "\"vbasis_" << i << "_y\" ";
	}
	outp << endl;

	while(y < 1) {
		double x = 0;
		while (x <= 1) {
			auto el_it = find_element(point(x, y, 0));
			if (el_it != nullptr) {
				
				outp << x << " " << y << " ";

				for (auto& dof_it : dofs) {
					vec3d res1(0, 0, 0);
					auto el_dof = el_it->get_dofs_num();
					for (size_t dof_i = 0; dof_i < el_dof.size(); dof_i++) {
							res1 = res1 + solutions[dof_it.number][el_dof[dof_i]] *  el_it->get_vector_basis_dof(dof_i)(x, y, 0);
					}
					outp << res1.x << " " << res1.y << " ";
				}

				outp << endl;
			}
			x += h;
		}
		y += h;
	}

	outp.close();

}

void VirtualFace_Hrot::set_mesh_files(string main_mesh, string bound_mesh) {
	main_mesh_file_name = main_mesh;
	bound_mesh_file_name = bound_mesh;
}


vector<dof_type> VirtualFace_Hrot::get_dofs_num() {
	vector<dof_type> dof_global_nums;
	dof_global_nums.reserve(dofs_n);
	for (auto& dof_i : dofs) {
		dof_global_nums.push_back(loc_to_glob[dof_i.number]);
	}

	return dof_global_nums;
}


VirtualFace_Hrot::~VirtualFace_Hrot() {
	for(auto& edge : edges)
		delete edge;
}

