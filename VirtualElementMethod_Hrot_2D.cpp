#include "VitualElementMethod_Hrot_2D.h"

VirtualElementMethod_Hrot_2D::VirtualElementMethod_Hrot_2D() {
	dofs_n = 1;
}

double VirtualElementMethod_Hrot_2D::get_lambda(VirtualFace_Hrot* el) {
	return 1.0;
}

double VirtualElementMethod_Hrot_2D::bound_func(dof_type basis_i, dof_type cur_row) {
	vec3d tau = local_edges[cur_row]->get_direction();
	double res = vec3d(0.0, 1.0, 0.0) * tau;

	return res;
}

void VirtualElementMethod_Hrot_2D::calculate() {
	// Расчитаем базичные функции на элементах
#pragma omp parallel for
	for (int el_i = 0; el_i < elements_n; el_i++) {
		elements[el_i]->input_mesh_from_params();
		elements[el_i]->calculate();
		#ifdef DEBUGOUTP
			cout << "Printing basis of element " << el_i << endl;
			stringstream db_ss;
			string db_basis_file_name;
			db_ss << "local_virtual_basis_of_" << el_i << ".txt";
			db_ss >> db_basis_file_name;
			elements[el_i]->test_print_local_basis(db_basis_file_name);
		#endif // DEBUGOUTP

		
	}
	vector<vfunc3d> rp_functions;
	rp_functions.push_back([&](double x, double y, double z)->vec3d {
		return get_lambda(nullptr)*vec3d(0, 1, 0);

	});


	generate_port();
	generate_matrix_with_out_bound(rp_functions);
	generate_matrix_first_bound([&](dof_type basis_i, dof_type cur_row)->double {
		return bound_func(basis_i, cur_row);
	});

	solve_SLAE();
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
	return nullptr;
}
