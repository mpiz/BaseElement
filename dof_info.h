#pragma once
#include <vector>
using namespace std;

typedef unsigned int dof_type;

namespace counters {
	static dof_type get_next_dof() {
		static dof_type dof_count = 0;
		return dof_count++;
	}
};

class dof_info {
 public:
	dof_info();
	dof_info(dof_type dof_n, dof_type ord = 0, dof_type nm = 0, dof_type cnt = 0);
	void set_geom(vector<dof_type>& geom_s);
	void add_geom(dof_type g);

	bool operator < (const dof_info& b) const;
	bool operator < (const dof_type& b) const;

	~dof_info();

	dof_type number;
	vector<dof_type> geom;

	dof_type order;
	dof_type num;
	dof_type count;

};

