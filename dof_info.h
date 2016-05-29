#pragma once
#include <vector>
#include <map>
using namespace std;

typedef unsigned int dof_type;

namespace counters {
	
	static dof_type get_next_global_dof() {
		static dof_type dof_count = 0;
		return dof_count++;
	}

	static dof_type get_next_dof(void* obj) {

		static map<void*, dof_type> local_counter;

		if (local_counter.find(obj) != local_counter.end()) {
			return local_counter[obj]++;
		}	else {
			local_counter[obj] = 1;
			return 0;
		}
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

