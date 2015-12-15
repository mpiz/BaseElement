#include "VirtualEdge_Hrot.h"


VirtualEdge_Hrot::VirtualEdge_Hrot() {
}


void VirtualEdge_Hrot::input_mesh(string file_name) {
	BaseElement::input_mesh(file_name, 102);
}

void VirtualEdge_Hrot::get_elements(vector<sector>& elems_g, int& dofs_n_g) {
	elems_g = elements;
	dofs_n_g = dofs_n;
}

void VirtualEdge_Hrot::set_right_parts(const vector<vbfunc>& right_parts_s) {
	right_parts = right_parts_s;
	dofs_n = right_parts.size();
}

void VirtualEdge_Hrot::calculate() {

	generate_matrix_with_out_bound(right_parts);
}

VirtualEdge_Hrot::~VirtualEdge_Hrot() {
}
