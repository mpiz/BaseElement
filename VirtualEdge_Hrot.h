#pragma once
#include "baseelement.h"
class VirtualEdge_Hrot : public BaseElement<sector> {
public:
	VirtualEdge_Hrot();

	void input_mesh(string file_name);
	void get_elements(vector<sector>& elems_g, int& dofs_n_g);
	void calculate();
	void set_right_parts(const vector<vfunc3d>& right_parts_s);

	~VirtualEdge_Hrot();

private:
	vector<vfunc3d> right_parts;
};

