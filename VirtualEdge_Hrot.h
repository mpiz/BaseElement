#pragma once
#include "baseelement.h"
class VirtualEdge_Hrot : public BaseElement<sector> {
public:
	VirtualEdge_Hrot();
	VirtualEdge_Hrot(dof_type order, dof_type num = 0);

	void input_mesh(string file_name);
	void get_elements(vector<sector*>& elems_g, int& dofs_n_g);
	void calculate();
	void set_right_parts(const vector<vfunc3d>& right_parts_s);

	vector<dof_info> calc_element_dofs(vector<node>& el_nodes);

	~VirtualEdge_Hrot();


	void test_calc_points(dof_type dof_i);

	dof_type get_dof_num(tuple<int, int, int> tup);
	vector<tuple<int, int, int>> get_bound_funcs();

private:
	vector<vfunc3d> right_parts;
};

