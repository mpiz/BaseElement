#pragma once

#include "BaseElement.h"


class VirtualFace_Hrot :	public BaseElement<trelement> {
public:
	VirtualFace_Hrot();

	VirtualFace_Hrot(const vector<node>& nodes_s);


	~VirtualFace_Hrot();

private:

	size_t edges_n;	// Количество ребёр

};

