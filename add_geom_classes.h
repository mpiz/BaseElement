#include "geom_classes.h"


class plane {
 public:
	 plane();
	 plane(const array<node, 3>& points);

	point to_local_cord(point p_glob);
	point to_global_cord(point p_loc);
	vec3d to_global_cord(vec3d v_loc);


	double get_jacobian() const;
	vec3d get_base_vec(size_t i) const;
	vec3d get_tau(size_t i) const;

	bool is_on_plane(point pn);

	static bool is_on_line(const array<node, 3>& points);

 private:

	 void init_cords(); //построение локальных координат
	 array<node, 3> plane_points; // Три точки, определяющие плоскость

	point trpoint[3]; //локальные координаты точек треугольника
	matrix(3) transition_matrix; //матрица перехода в локальные координаты

	array<vec3d, 4> tau;

	vec3d normal_vector;

	double jacobian;
	double h[2];
	array<vec3d, 2> base_vec;


};