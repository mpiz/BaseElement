#include "dof_info.h"

dof_info::dof_info() {
}

dof_info::dof_info(dof_type dof_n, dof_type ord, dof_type nm, dof_type cnt) {
	number = dof_n;
	order = ord;
	num = nm;
	count = cnt;

}

void dof_info::set_geom(vector<dof_type>& geom_s) {
	geom = geom_s;
}

void dof_info::add_geom(dof_type g) {
	geom.push_back(g);
}

bool dof_info::operator < (const dof_info& b) const {
	return number < b.number;
}

bool dof_info::operator < (const dof_type& b) const {
	return number < b;
}

dof_info::~dof_info() {
}
