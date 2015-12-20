#pragma once

#include "elements_classes.h"
#include "CGM.h"
#include <map>
#include <set>
#include <iostream>



template<typename elementT> class BaseElement {
 public:
	 BaseElement();
	 BaseElement(vector<node>& nodes_s, vector<dof_type>& dofs_s);

	 vector<dof_type> get_dofs();
	 vector<dof_type> get_local_dofs(); //Возвращает локальные степени свободы в глобальной нумерации

	virtual void calculate() = 0;
	virtual void input_mesh(string file_name) = 0;
	virtual vector<dof_info> calc_element_dofs(vector<node>& el_nodes); // Функци получения степеней свободы элемента по его узлам
	vector<dof_info> calc_element_dofs_edge(vector<node>& el_nodes, int dof_per_element); // Функци получения степеней свободы для рёбер

	virtual double vales(dof_type glob_dof_n, node pn);	// По глобальному номеру базисной функции и точке вычислем значение
	virtual double scalar_basis_v(dof_type loc_dof_n, node pn);	// По глобальному номеру базисной функции и точке вычислем значение

	virtual elementT* find_element(point pn);

	virtual double get_lambda(elementT* el);

	void input_mesh(string file_name, int valide_code); // ввод сетки из файла для виртуального элемента

	void generate_port();
	template<typename func_t> void generate_matrix_with_out_bound(vector<func_t> equation_right_part);
	void generate_matrix_first_bound(function<double(dof_type, dof_type)> bound_right_part); // Генерировать краевые условия, bound_right_part вектор вычисления краевых условий, 
																							 // Параметр 1 - номер функции, для которой считаем краевое
																							 // Параметр 2 - номер степени свободы, для которой считем краевое

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

	void print_full_matrix(string file_name);
	void print_right_part(dof_type dof_i, string file_name);

 protected:

	  void add_port(int add_el1, int add_el2);
	  int find_pos(int i, int j);

	 vector<elementT*> elements;

	 vector<vector<int>> port_colector;

	 vector<int> gi, gj;
	 vector<double> gg, di;

	 vector<double*> solutions;
	 vector<double*> rp;

	 vector<double> virtual_solution;

	 int elements_n;
	 size_t dofs_n;
	 size_t nodes_n;

	 size_t local_dof_n;
	 size_t local_nodes_n;

	 vector<dof_info> local_dofs;	// Степени свободы внутри элемента
	 vector<node> local_nodes;	// Узлы внутри элемента

	 vector<dof_info> dofs;		// Степени свободы для всей виртуальной сетки
	 vector<node> nodes;	// Узлы для всей виртуальной сетки


	 map<dof_type, dof_type> virtual_glob_loc;
	 map<dof_type, dof_type> loc_to_glob;
	 map<dof_type, dof_type> glob_to_loc;

	 vector<dof_type> first_bound;

	 vector<func3d> basis_right_parts;

	 CGM solver;

	 map<tuple<int, int, int>, dof_info> edge_type_dofs; // Edge-степени свободы для векторных методов первые два int - номера узлов, второй - номер функции
	 map<dof_info, tuple<int, int, int>> edge_type_dofs_revers; // Реверсированный массив edge_type_dofs
	 dof_type local_dof_counter;

	 dof_type method_order, method_num;
};

// ========================================================================
// =========================== Реализации =================================
// ========================================================================


template<typename elementT> BaseElement<elementT>::BaseElement() {
	local_dof_counter = 0;
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

template<typename elementT> double BaseElement<elementT>::get_lambda(elementT* el) {
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
	int cur = -1;
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
		vector<dof_type> loc_dof = elements[el_i]->get_dofs_num();
		int loc_dof_n = loc_dof.size();
		
		for(int i = 0; i < loc_dof_n; i++)
			for(int j = 0; j < i; j++)
				add_port(loc_dof[i], loc_dof[j]);
	}

	//Собираем портрет
	int gi_it = 0;
	gi[gi_it] = 0;
	gi_it++;


	for(size_t port_i = 0; port_i < local_dof_n; port_i++) {
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

	for(size_t i = 0; i < dofs_n; i++) {
		rp[i] = new double [local_dof_n];
		solutions[i] = new double [local_dof_n];
	}

	// Обнуление
	
	for(size_t i = 0; i < local_dof_n; i++) {
		di[i] = 0;
		for(size_t k = 0; k < dofs_n; k++) {
			rp[k][i] = 0;
			solutions[k][i] = 0;
		}
	}

	for(size_t i = 0; i < gg.size(); i++) 
		gg[i] = 0;

	// Собрка
	for(int el_i = 0; el_i < elements_n; el_i++) {
		auto el_dof = elements[el_i]->get_dofs();
		size_t el_dof_n = el_dof.size();

		double lambda = get_lambda(elements[el_i]);

		auto A_loc = elements[el_i]->get_local_matrix(lambda);

		// Для каждого уравнения будет своя правая часть
		for(size_t k = 0; k < dofs_n; k++) {
			auto b_loc = elements[el_i]->get_local_right_part(equation_right_part[k]);
			for(size_t i = 0; i < el_dof_n; i++)
				rp[k][el_dof[i].number] += b_loc[i];
		}
		

		for(size_t i = 0; i < el_dof_n; i++) {
			int i_dof = el_dof[i].number;
			for(size_t j = 0; j < i; j++) {
				int pos = find_pos(i_dof,el_dof[j].number);
				gg[pos] += A_loc[i][j];
			}
			di[i_dof] += A_loc[i][i];

		}

	}


}

template<typename elementT> void BaseElement<elementT>::generate_matrix_first_bound(function<double(dof_type, dof_type)> bound_right_part) {
	vector<double> vals;
	vals.resize(dofs_n);

	auto fb_size = first_bound.size();

	for(size_t k = 0; k < fb_size; k++)	{
		int cur_row = first_bound[k];

		for(size_t basis_i = 0; basis_i < dofs_n; basis_i++) {

			vals[basis_i] = bound_right_part(basis_i, cur_row);
			rp[basis_i][cur_row] = vals[basis_i];
		}


		di[cur_row] = 1;

		int i_s = gi[cur_row], i_e = gi[cur_row+1];
		for(int i = i_s; i < i_e; i++){
			for(size_t basis_i = 0; basis_i < dofs_n; basis_i++) {
				rp[basis_i][gj[i]] -= gg[i]*vals[basis_i];
			}
			gg[i] = 0;
		}
		for(size_t p = cur_row + 1; p < local_dof_n; p++){
			int i_s = gi[p], i_e = gi[p+1];
			for(int i = i_s; i < i_e; i++){
				if(gj[i] == cur_row){
					for(size_t basis_i = 0; basis_i < dofs_n; basis_i++)
						rp[basis_i][p] -= gg[i]*vals[basis_i];
					gg[i] = 0;
				}
			}
		}
	}


}

template<typename elementT> void BaseElement<elementT>::solve_SLAE() {
	for(size_t basis_i = 0; basis_i < dofs_n; basis_i++) {
			solver.init(gi.data(), gj.data(), di.data(), gg.data(), local_dof_n);
			solver.solve(rp[basis_i], solutions[basis_i]);

#ifdef DEBUGOUTP
	auto sol = solutions[basis_i];
	continue;
#endif 
	}
}

template<typename elementT> void BaseElement<elementT>::print_full_matrix(string file_name) {

	ofstream outp(file_name.c_str());

	for(size_t i = 0; i < local_dof_n; i++) {
		for(size_t j = 0; j < local_dof_n; j++) {
			double val;
			if (i == j) 
				val = di[i];
			else {
				size_t k = find_pos(i, j);
				if( k == -1)
					val = 0;
				else
					val = gg[k];
			}

				outp << val;
				if (j != local_dof_n-1)
					outp << "\t";
		}
		outp << endl;
	}

	outp.close();

}

template<typename elementT> void BaseElement<elementT>::print_right_part(dof_type dof_i, string file_name) {

	ofstream outp(file_name.c_str());

	for(size_t i = 0; i < local_dof_n; i++) {
		outp << rp[dof_i][i] << endl;
	}

	outp.close();

}

template<typename elementT> vector<dof_info> BaseElement<elementT>::calc_element_dofs_edge(vector<node>& el_nodes, int dof_per_element) {
	vector<dof_info> el_dofs;
	auto el_nodes_n = el_nodes.size();

	set<dof_type> tmp_dofs; // Извращение нужно для отрезков, т.к. у них на 2 узла одно ребро, у треугольников и тетраэдров с этим всё ок

	for(size_t node_i = 0; node_i < el_nodes_n; node_i++) {
		size_t next_node = (node_i + 1) % el_nodes_n;

		int n1 = std::min(el_nodes[node_i].number, el_nodes[next_node].number);
		int n2 = std::max(el_nodes[node_i].number, el_nodes[next_node].number);
		for(int k = 0; k < dof_per_element; k++) {
			auto tup = make_tuple(n1, n2, k);
			auto map_res = edge_type_dofs.find(tup);
			dof_info ins_dof;
			if (map_res != edge_type_dofs.end()) {
				ins_dof = map_res->second;
			}
			else {
				ins_dof = dof_info(local_dof_counter, 1, k);
				ins_dof.add_geom(n1);
				ins_dof.add_geom(n2);
				local_dof_counter++;

				local_dof_n = local_dof_counter;
				local_dofs.push_back(ins_dof);
				edge_type_dofs[tup] = ins_dof;
				edge_type_dofs_revers[ins_dof] = tup;
			}

			if(tmp_dofs.find(ins_dof.number) == tmp_dofs.end()) {
				el_dofs.push_back(ins_dof);
				tmp_dofs.insert(ins_dof.number);
			}

		}
	}

	return el_dofs;
}

template<typename elementT> vector<dof_info> BaseElement<elementT>::calc_element_dofs(vector<node>& el_nodes) {
	// В простейшем случае степени свободы - номера узлов
	vector<dof_info> el_dofs;
	for(auto& node_it : el_nodes) {
		dof_info d(node_it.number);
		d.add_geom(node_it.number);
		el_dofs.push_back(d);
	}

	return el_dofs;
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

	local_dof_n = local_nodes_n;

	vector<dof_info> tr_dofs;
	vector<node> tr_node;

	// Вводим узлы
	for(size_t i = 0; i < local_nodes_n; i++) {
		double x, y, z;
		int tmp_int;

		inp_file >> tmp_int >> x >> y >> z;
		tmp_int--;
		loc_to_glob[i] = tmp_int;
		glob_to_loc[tmp_int] = i;
		local_nodes[i] = node(x,y,z);
		local_nodes[i].number = tmp_int;
	}

	int element_nodes = elementT::element_nodes;

	//Вводим элементы
	for(int i = 0; i < total_element_n; i++) {
		int num, el_code, tr_p;

		tr_dofs.clear();
		tr_node.clear();

		inp_file >> num >> el_code;
		for(int j = 0; j < dat_codes[el_code]; j++) {
			inp_file >> tr_p;
			tr_p--;
			tr_p = glob_to_loc[tr_p];
			tr_node.push_back(local_nodes[tr_p]);
		}
		tr_dofs = calc_element_dofs(tr_node);
		if(el_code == valide_code)
			elements.push_back(new elementT(tr_node, tr_dofs));
	
	}
	elements_n = elements.size();

}