#include "VirtualEdge_Hrot.h"


VirtualEdge_Hrot::VirtualEdge_Hrot() {
}

VirtualEdge_Hrot::VirtualEdge_Hrot(dof_type order, dof_type num) {
	method_order = order;
	method_num = num;
}


void VirtualEdge_Hrot::input_mesh(string file_name) {
	BaseElement::input_mesh(file_name, 102);
}

void VirtualEdge_Hrot::get_elements(vector<sector>& elems_g, int& dofs_n_g) {
	elems_g = elements;
	dofs_n_g = dofs_n;
}

void VirtualEdge_Hrot::set_right_parts(const vector<vfunc3d>& right_parts_s) {
	right_parts = right_parts_s;
	dofs_n = right_parts.size();
}

vector<dof_type> VirtualEdge_Hrot::calc_element_dofs(vector<node>& el_nodes) {
	return calc_element_dofs_edge(el_nodes, sector::get_dof_n(method_order, method_num));
}

void VirtualEdge_Hrot::calculate() {

	generate_port();
	generate_matrix_with_out_bound(right_parts);
	solve_SLAE();
}

void VirtualEdge_Hrot::test_calc_points(dof_type dof_i) {
	struct test_p {
		double x, y, z;
		double val;
		vec3d vval;
	};

	vector<test_p> vals;

	for(auto& sect : elements) {
		auto loc_dof = sect.get_dofs();
		sect.for_point_on_element([&](double x, double y, double z){
			test_p tp;
			tp.x = x; tp.y = y; tp.z = z;
			tp.val = 0;
			tp.vval = vec3d(0, 0, 0);
			for(int i = 0; i < loc_dof.size(); i++) {
				double q = solutions[dof_i][loc_dof[i]] ;
				double tau = sect.get_vector_basis_dof_tau(i)(x, y, z);
				vec3d v = sect.get_vector_basis_dof(i)(x, y, z);//right_parts[dof_i](x, y, z);//
				tp.val += q * tau; 

				tp.vval = tp.vval + q * v;
			}

			vals.push_back(tp);
		});

	}

	ofstream outp("test_1.txt");
	for(auto& t : vals) {
		outp << t.x << " " << t.y << " " << t.vval.x << " " << t.vval.y << endl;
	}
	outp.close();

}

VirtualEdge_Hrot::~VirtualEdge_Hrot() {
}
