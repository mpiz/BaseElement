#include "add_geom_classes.h"

plane::plane() {

}

plane::plane(const array<node, 3>& points) {
	plane_points = points;
	init_cords();
}

double plane::get_jacobian() const {
	return jacobian;
}

vec3d plane::get_base_vec(size_t i) const {
	return base_vec[i];
}

vec3d plane::get_tau(size_t i) const {
	return tau[i];
}

point plane::to_local_cord(point p_glob) {
	point p_shift = p_glob;

	//сдвиг

	for(int i = 0; i < 3; i++)
		p_shift[i] -= plane_points[0][i];

	point p_loc;

	//поворот
	for(int i = 0; i < 3; i++) {
		p_loc[i] = 0;
		for(int j = 0; j < 3; j++)
			p_loc[i] += transition_matrix[i][j] * p_shift[j];
	}

	return p_loc;
}

point plane::to_global_cord(point p_loc) {
	point p_turn;

	//поворот
	for(int i = 0; i < 3; i++) {
		p_turn[i] = 0;
		for(int j = 0; j < 3; j++)
			p_turn[i] += transition_matrix[j][i] * p_loc[j];
	}

	point p_glob; 

	//сдвиг
	for(int i = 0; i < 3; i++)
		p_glob[i] = p_turn[i] + plane_points[0][i];

	return p_glob;

}

vec3d plane::to_global_cord(vec3d v_loc) {
	vec3d v_glob;

	//поворот
	for(int i = 0; i < 3; i++) {
		v_glob[i] = 0;
		for(int j = 0; j < 3; j++)
			v_glob[i] += transition_matrix[j][i] * v_loc[j];
	}


	return v_glob;

}

bool plane::is_on_line(const array<node, 3>& points) {

	vec3d line_a = vec3d(points[0], points[1]);
	vec3d line_b = vec3d(points[0].x, points[0].y, points[0].z);

	double length = line_a.norm();
	double length_sq = length * length;
	size_t line_num;

	if (fabs(line_a.x) > 1e-16)
		line_num = 0;
	else
		if (fabs(line_a.y) > 1e-16)
			line_num = 1;
		else
			line_num = 2;

	double t = (points[2][line_num] - line_b[line_num]) / line_a[line_num];
	bool cond1 = (t >= -geom_err && t <= (1 + geom_err));
	bool line_cond[3];
	bool cond2 = true;
	for (int i = 0; i < 3; i++) {
		double s = line_a[i] * t + line_b[i] - points[2][i];
		line_cond[i] = (fabs(line_a[i] * t + line_b[i] - points[2][i]) < geom_err);
		cond2 = cond2 && line_cond[i];
	}
	return cond1 && cond2;
}


bool plane::is_on_plane(point pn) {
	return true;
};

void plane::init_cords() {
		//Построение локальной системы координат
	base_vec[0] = vec3d(plane_points[0], plane_points[1]); 
	base_vec[1] = vec3d(plane_points[0], plane_points[2]); 
	vec3d e2 = base_vec[1] - ((base_vec[0]*base_vec[1]) / (base_vec[0]*base_vec[0])) * base_vec[0];
	vec3d e3 = base_vec[0].cross(e2);
	normal_vector = e3;

	tau[0] = base_vec[0];
	tau[1] = base_vec[1];
	tau[2] = vec3d(plane_points[1], plane_points[2]);

	h[0] = base_vec[0].norm();
	h[1] = base_vec[1].norm();

	trpoint[0] = point(0, 0, 0);
	trpoint[1] = point(base_vec[0].norm(), 0, 0);

	vec3d e1 = base_vec[0] / base_vec[0].norm();
	e2 = e2 / e2.norm();
	e3 = e3 / e3.norm();

	for(int i = 0; i < 3; i++) {
		transition_matrix[0][i] = e1[i];
		transition_matrix[1][i] = e2[i];
		transition_matrix[2][i] = e3[i];
	}

	trpoint[2] = (transition_matrix * base_vec[1]).to_point();
	jacobian = (base_vec[0].cross(base_vec[1])).norm();
}