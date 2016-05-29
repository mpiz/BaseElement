#include "VitualElementMethod_Hrot_2D.h"

//#define DEBUGOUTP_1

vec3d virtual_solution_func(double x, double y, double z) {
	return vec3d(0, 1, 0);
}

VirtualElementMethod_Hrot_2D::VirtualElementMethod_Hrot_2D() {
	dofs.push_back(dof_info(0));
	dofs_n = 1;
}

double VirtualElementMethod_Hrot_2D::get_lambda(VirtualFace_Hrot* el) {
	return 1.0;
}

double VirtualElementMethod_Hrot_2D::bound_func(dof_type basis_i, dof_type cur_row) {
	vec3d tau = local_edges[cur_row]->get_direction();

	auto edge_nodes = local_edges[cur_row]->get_nodes();
	auto edge_nodes_n = edge_nodes.size();
	point pn(0, 0, 0);
	for (int i = 0; i < edge_nodes_n; i++)
		for (int j = 0; j < 3; j++)
			pn[j] += edge_nodes[i][j] / edge_nodes_n;

	double res = virtual_solution_func(pn.x, pn.y, pn.z) * tau;

	return res;
}

void VirtualElementMethod_Hrot_2D::calculate() {
	// Расчитаем базичные функции на элементах

#pragma omp parallel for
	for (int el_i = 0; el_i < elements_n; el_i++) {
		cout << "Calulating " << el_i << endl;
		elements[el_i]->input_mesh_from_params();
		elements[el_i]->calculate();
		cout << "Calulating " << el_i << " finished" <<endl;
		#ifdef DEBUGOUTP_1
			cout << "Printing basis of element " << el_i << endl;
			stringstream db_ss;
			string db_basis_file_name;
			db_ss << "local_virtual_basis_of_" << el_i << ".txt";
			db_ss >> db_basis_file_name;
			elements[el_i]->test_print_local_basis(db_basis_file_name);
		#endif // DEBUGOUTP_1

		
	}
	vector<vfunc3d> rp_functions;
	rp_functions.push_back([&](double x, double y, double z)->vec3d {
		return get_lambda(nullptr)*virtual_solution_func(x, y, z);

	});

	cout << "Assembling matrix\n";
	generate_port();
	generate_matrix_with_out_bound(rp_functions);
	cout << "Assembling bounds\n";
	generate_matrix_first_bound([&](dof_type basis_i, dof_type cur_row)->double {
		return bound_func(basis_i, cur_row);
	});

#ifdef DEBUGOUTP_1
	print_full_matrix("virtual_matrix.txt");
	print_right_part(0, "vurtial_rp.txt");
#endif

	cout << "Solving SLAE\n";
	solve_SLAE();

#ifdef DEBUGOUTP
	FILE* foutp = fopen("SLAE_solution.txt", "w");
	for (int i = 0; i < local_dof_n; i++)
		fprintf(foutp, "%.15lf\n", solutions[0][i]);
	fclose(foutp);
#endif
}

void VirtualElementMethod_Hrot_2D::input_mesh(string file_name) {
	ifstream inp_file(file_name.c_str());

	// Введём количество узлов, рёбер, элементов и степеней свободы
	inp_file >> local_nodes_n >> local_edges_n >> elements_n;

	local_dof_n = local_edges_n;

	local_nodes.reserve(local_nodes_n);
	local_edges.reserve(local_edges_n);
	elements.reserve(local_edges_n);

	// Введём узлы
	for (int n_i = 0; n_i < local_nodes_n; n_i++) {
		double x, y;
		inp_file >> x >> y;
		local_nodes.push_back(node(x, y, 0, n_i));

	}

	// Введём рёбра
	for (int e_i = 0; e_i < local_edges_n; e_i++) {
		int n1, n2;
		vector<node> edge_nodes;
		edge_nodes.resize(2);
		inp_file >> n1 >> n2;
		if (n1 > n2)
			swap(n1, n2);

		edge_nodes[0] = local_nodes[n1];
		edge_nodes[1] = local_nodes[n2];


		dof_info edge_dof(e_i);
		edge_dof.add_geom(n1);
		edge_dof.add_geom(n2);
		vector<dof_info> edge_dofs;
		edge_dofs.push_back(edge_dof);
		
		nodes_to_edges[pair<size_t, size_t>(n1, n2)] = e_i;

		local_edges.push_back(new sector(edge_nodes, edge_dofs));
	}

	// Соберём элемент
	for (int el_i = 0; el_i < elements_n; el_i++) {
		int el_nodes_n;
		vector<node> el_nodes;

		inp_file >> el_nodes_n;
		el_nodes.reserve(el_nodes_n);
		for (int el_n_i = 0; el_n_i < el_nodes_n; el_n_i++) {
			int n_val;
			inp_file >> n_val;
			el_nodes.push_back(local_nodes[n_val]);
		}
		string main_mesh, bound_mesh;
		inp_file >> main_mesh >> bound_mesh;

		VirtualFace_Hrot* el = new VirtualFace_Hrot(el_nodes, nodes_to_edges, virtual_method_order);
		el->set_mesh_files(main_mesh, bound_mesh);
		elements.push_back(el);
	}



}

void VirtualElementMethod_Hrot_2D::input_bound(string file_name) {
	ifstream inp_file(file_name.c_str());
	int fb_size;
	inp_file >> fb_size;
	first_bound.reserve(fb_size);
	for (int i = 0; i < fb_size; i++) {
		int tmp;
		inp_file >> tmp;
		first_bound.push_back(tmp);

	}
}

vector<dof_info> VirtualElementMethod_Hrot_2D::calc_element_dofs(vector<node>& el_nodes) {
	return vector<dof_info>();
}

vec3d VirtualElementMethod_Hrot_2D::vector_basis_val(dof_type basis_i, double x, double y, double z) {
	return vec3d();
}

VirtualFace_Hrot * VirtualElementMethod_Hrot_2D::find_element(point pn) {
	for (auto el_i : elements) {
		if (el_i->find_element(pn) != nullptr)
			return el_i;
	}

	return nullptr;
}


void VirtualElementMethod_Hrot_2D::test_print_solution(string file_name) {
	ofstream outp(file_name.c_str());

	double y = 0;
	double h = 0.02;
	outp << "VARIABLES = \"x\" \"y\" ";
	for (auto i = 0; i < dofs_n; i++) {
		outp << "\"vbasis_" << i << "_x\" ";
		outp << "\"vbasis_" << i << "_y\" ";
		outp << "\"truesol_" << i << "_x\" ";
		outp << "\"truesol_" << i << "_y\" ";
		outp << "\"diff_" << i << "_x\" ";
		outp << "\"diff_" << i << "_y\" ";
		outp << "\"tau-11_" << i << "\" ";
	}
	outp << endl;
	
	vec3d tau_11(-1.0 / sqrt(2.0), 1 / sqrt(2.0), 0);

	while (y < 1) {
		double x = 0;
		while (x <= 1) {
			auto el_it = find_element(point(x, y, 0));
			if (el_it != nullptr) {

				outp << x << " " << y << " ";

				for (auto& dof_it : dofs) {
					vec3d res1(0, 0, 0);
					auto el_dof = el_it->get_dofs_num();
					for (dof_type dof_i = 0; dof_i < el_dof.size(); dof_i++) {
						res1 = res1 + solutions[dof_it.number][el_dof[dof_i]] * el_it->get_vector_basis_dof(dof_i, x, y, 0);
					}
					vec3d sol = virtual_solution_func(x, y, 0);
					vec3d diff = sol - res1;
					outp << res1.x << " " << res1.y << " " << sol.x << " " << sol.y << " " << diff.x << " " << diff.y << " " << res1*tau_11 << " ";
				}

				outp << endl;
			}
			x += h;
		}
		y += h;
	}

	outp.close();
}

vec3d VirtualElementMethod_Hrot_2D::func_in_element(size_t el_i, double x, double y, double z) {
	vec3d res(0, 0, 0);
	auto el_dofs = elements[el_i]->get_dofs_num();
	auto el_dofs_n = el_dofs.size();
	for (dof_type i = 0; i < el_dofs_n; i++) {
		res = res + solutions[0][el_dofs[i]] * elements[el_i]->get_vector_basis_dof(i, x, y, z);
	}

	return res;
}

double VirtualElementMethod_Hrot_2D::calc_error() {
	double res = 0;
	double f_val = 0;
	for (int el_i = 0; el_i < elements_n; el_i++) {
		res += elements[el_i]->integrate([&](double x, double y, double z)->double {
			vec3d sol_val = func_in_element(el_i, x, y, z);
			vec3d func_val = virtual_solution_func(x, y, z);
			vec3d diff = sol_val - func_val;
			double int_res = diff*diff;
			return int_res;

		});

		f_val += elements[el_i]->integrate([&](double x, double y, double z)->double {
			vec3d func_val = virtual_solution_func(x, y, z);
			return func_val * func_val;
		});
	}
	res = sqrt(res / f_val);
	return res;
}