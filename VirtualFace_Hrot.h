#pragma once

#include "BaseElement.h"


class VirtualFace_Hrot : public BaseElement<trelement> {
public:
	VirtualFace_Hrot();

	VirtualFace_Hrot(const vector<node>& nodes_s, const vector<dof_type>& dofs_s, dof_type order, dof_type num = 0);


	~VirtualFace_Hrot();

private:

	plane face_plane;
	size_t edges_n;	// Количество ребёр

	vector<sector> edges;
	size_t edges_n;

	vector<vbfunc> bound_functions;	// Функции для учёта краевых условий




};

