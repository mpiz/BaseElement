#pragma once

#include "elements_classes.h"
#include "CGM.h"
#include <map>
#include <set>
#include <iostream>

namespace counters {
	static dof_type get_next_dof() {
		static dof_type dof_count = 0;
		return dof_count++;
	}
};


template<typename elementT> class BaseElement {
 public:
	 BaseElement();
	 BaseElement(vector<node>& nodes_s, vector<dof_type>& dofs_s);

	 vector<dof_type> get_dofs();
	 vector<dof_type> get_local_dofs(); //Возвращает локальные степени свободы в глобальной нумерации

	virtual void calculate() = 0;
	virtual void input_mesh(string file_name) = 0;

	virtual double vales(dof_type glob_dof_n, node pn);	// По глобальному номеру базисной функции и точке вычислем значение
	virtual double scalar_basis_v(dof_type loc_dof_n, node pn);	// По глобальному номеру базисной функции и точке вычислем значение

	virtual elementT* find_element(point pn);

	virtual double get_lambda(elementT& el);

	void input_mesh(string file_name, int valide_code); // ввод сетки из файла для виртуального элемента

	void generate_port();
	template<typename func_t> void generate_matrix_with_out_bound(vector<func_t> equation_right_part);

	void solve_SLAE();

	void get_solutions(vector<double*>& solutions_g);
	void get_local_glob(map<dof_type, dof_type>& loc_to_glob_g);
	void get_glob_local(map<dof_type, dof_type>& glob_to_loc_g);
	void get_virtual_glob_local(map<dof_type, dof_type>& virtual_glob_to_loc_g);
	void get_basis_right_parts(vector<func3d>& basis_right_parts_g);

	double value_element(elementT* elem, point pn, dof_type dof_i = 0);
	
	int get_local_elements_n();
	int get_local_nodes_n();
	int get_local_dofs_n();

	int& operator [] (int loc_n);

	void set_virtual_solution(vector<double>& virt_sol);
	void set_virtual_solution_to_elements();

 protected:

	  void add_port(int add_el1, int add_el2);
	  int find_pos(int i, int j);

	 vector<elementT> elements;

	 vector<vector<int>> port_colector;

	 vector<int> gi, gj;
	 vector<double> gg, di;

	 vector<double*> solutions;
	 vector<double*> rp;

	 vector<double> virtual_solution;

	 int elements_n;
	 int dofs_n;
	 int nodes_n;

	 int local_dof_n;
	 int local_nodes_n;

	 vector<dof_type> local_dofs;	// Степени свободы внутри элемента
	 vector<node> local_nodes;	// Узлы внутри элемента

	 vector<dof_type> dofs;		// Степени свободы для всей виртуальной сетки
	 vector<node> nodes;	// Узлы для всей виртуальной сетки


	 map<dof_type, dof_type> virtual_glob_loc;
	 map<dof_type, dof_type> loc_to_glob;
	 map<dof_type, dof_type> glob_to_loc;

	 vector<dof_type> first_bound;

	 vector<func3d> basis_right_parts;

	 CGM solver;

	 size_t	element_order;

};

// ========================================================================
// =========================== Реализации =================================
// ========================================================================


template<typename elementT> BaseElement<elementT>::BaseElement() {

}

template<typename elementT> BaseElement<elementT>::BaseElement(vector<node>& nodes_s, vector<dof_type>& dofs_s) {
	nodes = nodes_s;
	dofs = dofs_s;

	nodes_n = nodes.size();
	dofs_n = dofs.size();
	virtual_solution.resize(dofs_n);

	for(int i = 0; i < dofs_n; i++) {
		virtual_glob_loc[dofs[i]] = i;
	}
}

template<typename elementT> vector<dof_type> BaseElement<elementT>::get_dofs() {
	return dofs;
}
template<typename elementT> vector<dof_type> BaseElement<elementT>::get_local_dofs() {
	return local_dofs;
}

template<typename elementT> void BaseElement<elementT>::get_solutions(vector<double*>& solutions_g) {
	solutions_g = solutions;
}

template<typename elementT> void BaseElement<elementT>::get_local_glob(map<dof_type,dof_type>& loc_to_glob_g){
	loc_to_glob_g = loc_to_glob;
}

template<typename elementT> void BaseElement<elementT>::get_glob_local(map<dof_type,dof_type>& glob_to_loc_g){
	glob_to_loc_g = glob_to_loc;
}

template<typename elementT> void BaseElement<elementT>::get_virtual_glob_local(map<dof_type,dof_type>& virtual_glob_to_loc_g){
	virtual_glob_to_loc_g = virtual_glob_loc;
}

template<typename elementT> void BaseElement<elementT>::get_basis_right_parts(vector<func3d>& basis_right_parts_g){
	basis_right_parts_g = basis_right_parts;
}

template<typename elementT> int BaseElement<elementT>::get_local_elements_n() {
	return elements_n;
}

template<typename elementT> int BaseElement<elementT>::get_local_nodes_n() {
	return local_nodes_n;
}

template<typename elementT> int BaseElement<elementT>::get_local_dofs_n() {
	return local_dof_n;
}

template<typename elementT> int& BaseElement<elementT>::operator[] (int loc_n) {
	return dofs[loc_n];
}

template<typename elementT> double BaseElement<elementT>::get_lambda(elementT& el) {
	return 1.0;
}

template<typename elementT> void BaseElement<elementT>::set_virtual_solution(vector<double>& virt_sol) {
	int virt_sol_size = virt_sol.size();
	virtual_solution.clear();
	for(int i = 0; i < virt_sol_size; i++)
		virtual_solution.push_back(virt_sol[i]);
}

template<typename elementT> double BaseElement<elementT>::vales(dof_type loc_dof_n, node pn) {
	return 0;
}
template<typename elementT> double BaseElement<elementT>::scalar_basis_v(dof_type loc_dof_n, node pn) {
	return 0;
}

template<typename elementT> elementT* BaseElement<elementT>::find_element(point pn) {
	throw;
	return nullptr;
}

template<typename elementT> double BaseElement<elementT>::value_element(elementT* elem, point pn, dof_type dof_i = 0) {
	if (elem == nullptr)
		return 0;

	auto local_dofs = elem->get_dofs();
	int local_dofs_size = local_dofs.size();
	double value = 0;

#ifdef DEBUGOUTP
	vector<double> basis_vals, q_vals;
	basis_vals.resize(local_dofs_size);
	q_vals.resize(local_dofs_size);
#endif

	for (int i = 0; i < local_dofs_size; i++) {
		double basis_v = elem->scalar_basis_v(i, pn.x, pn.y, pn.z);
		double q_i = solutions[dof_i][local_dofs[i]];
		#ifdef DEBUGOUTP
		basis_vals[i] = basis_v;
		q_vals[i] = q_i;
		#endif
		value += q_i * basis_v;
	}

	return value;
	
}

template<typename elementT> int BaseElement<elementT>::find_pos(int i, int j) {
	if(j > i)
		swap(i, j);

	int k_s = gi[i], k_e = gi[i+1];
	int cur;
	bool find = false;
	for(int k = k_s; k < k_e && !find; k++){
		if(gj[k] == j){
			cur = k;
			find = true;
		}
	}
	return cur;
}


template<typename elementT> void BaseElement<elementT>::add_port(int add_el1, int add_el2) {

	
	int what, where; //where

	if(add_el1 > add_el2) {
		what = add_el2;
		where = add_el1;
	}
	else {
		what = add_el1;
		where = add_el2;
	}
	bool is_there = false;

	for(auto iter = port_colector[where].begin(); iter != port_colector[where].end() && !is_there; iter++) {
		if((*iter) == what)
			is_there = true;
	}

	if(!is_there){
		port_colector[where].push_back(what);
	}

}


template<typename elementT> void BaseElement<elementT>::generate_port() {

	port_colector.resize(local_dof_n);
	gi.resize(local_dof_n + 1);

	//Собираем портрет
	for(int el_i = 0; el_i < elements_n; el_i++) {
		vector<dof_type> loc_dof = elements[el_i].get_dofs();
		int loc_dof_n = loc_dof.size();
		
		for(int i = 0; i < loc_dof_n; i++)
			for(int j = 0; j < i; j++)
				add_port(loc_dof[i], loc_dof[j]);
	}

	//Собираем портрет
	int gi_it = 0;
	gi[gi_it] = 0;
	gi_it++;


	for(int port_i = 0; port_i < local_dof_n; port_i++) {
		sort(port_colector[port_i].begin(), port_colector[port_i].end());
		gi[gi_it] = gi[gi_it-1] + port_colector[port_i].size();
		gi_it++;
		for(auto iter1 = port_colector[port_i].begin(); iter1 != port_colector[port_i].end(); iter1++) {
			gj.push_back(*iter1);
		}
	}

	gg.resize(gj.size());
	di.resize(local_dof_n);
	port_colector.resize(0);
	
}



template<typename elementT> template<typename func_t> void BaseElement<elementT>::generate_matrix_with_out_bound(vector<func_t> equation_right_part) {
	rp.resize(dofs_n);
	solutions.resize(dofs_n);

	for(int i = 0; i < dofs_n; i++) {
		rp[i] = new double [local_dof_n];
		solutions[i] = new double [local_dof_n];
	}

	// Обнуление
	
	for(int i = 0; i < local_dof_n; i++) {
		di[i] = 0;
		for(int k = 0; k < dofs_n; k++) {
			rp[k][i] = 0;
			solutions[k][i] = 0;
		}
	}

	for(int i = 0; i < gg.size(); i++) 
		gg[i] = 0;

	// Собрка
	for(int el_i = 0; el_i < elements_n; el_i++) {
		auto el_dof = elements[el_i].get_dofs();
		int el_dof_n = el_dof.size();

		double lambda = get_lambda(elements[el_i]);

		auto A_loc = elements[el_i].get_local_matrix(lambda);

		// Для каждого уравнения будет своя правая часть
		for(int k = 0; k < dofs_n; k++) {
			auto b_loc = elements[el_i].get_local_right_part(equation_right_part[k]);
			for(int i = 0; i < el_dof_n; i++)
				rp[k][el_dof[i]] += b_loc[i];
		}
		

		for(int i = 0; i < el_dof_n; i++) {
			int i_dof = el_dof[i];
			for(int j = 0; j < i; j++) {
				int pos = find_pos(i_dof,el_dof[j]);
				gg[pos] += A_loc[i][j];
			}
			di[i_dof] += A_loc[i][i];

		}

		// Базис первого порядка - правая часть нулевая
	}


}

template<typename elementT> void BaseElement<elementT>::solve_SLAE() {
	for(int basis_i = 0; basis_i < dofs_n; basis_i++) {
			solver.init(&gi.front(), &gj.front(), &di.front(), &gg.front(), local_dof_n);
			solver.solve(rp[basis_i], solutions[basis_i]);
	}
}


template<typename elementT> void BaseElement<elementT>::input_mesh(string file_name, int valide_code) {
	ifstream inp_file(file_name.c_str());

	map<int, int> dat_codes;
	dat_codes[102] = 2; // отрезок
	dat_codes[203] = 3; // треугольник
	dat_codes[304] = 4; // тетраэдр
	int total_element_n;

	inp_file >> local_nodes_n >> total_element_n;
	local_nodes.resize(local_nodes_n);
	local_dofs.resize(local_nodes_n);

	local_dof_n = local_nodes_n;

	vector<dof_type> tr_dofs;
	vector<node> tr_node;

	// Вводим узлы
	for(int i = 0; i < local_nodes_n; i++) {
		double x, y, z;
		int tmp_int;

		inp_file >> tmp_int >> x >> y >> z;
		tmp_int--;
		loc_to_glob[i] = tmp_int;
		glob_to_loc[tmp_int] = i;
		local_nodes[i] = node(x,y,z);
		local_nodes[i].number = tmp_int;
		local_dofs[i] = tmp_int;
	}

	int element_nodes = elementT::element_nodes;

	//Вводим элементы
	for(int i = 0; i < total_element_n; i++) {
		int num, el_code, tr_p;
		int tmp_int;

		tr_dofs.clear();
		tr_node.clear();

		inp_file >> num >> el_code;
		for(int j = 0; j < dat_codes[el_code]; j++) {
			inp_file >> tr_p;
			tr_p--;
			tr_p = glob_to_loc[tr_p];
			tr_node.push_back(local_nodes[tr_p]);
			tr_dofs.push_back(tr_p);
		}
		if(el_code == valide_code)
			elements.push_back(elementT(tr_node, tr_dofs));
	
	}
	elements_n = elements.size();

}