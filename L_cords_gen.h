#pragma once
//���������� L-��������� ��� ������������ � ���������

#include "macro.h"

//========= ��������� ==========

double determenant3(matrix(3) A); //������������ ������� 3x3
matrix(4) inverse4(matrix(4) A, double& detA); //��������� ������� 4x4 (����� �������)

//========= ������������ ==========

matrix(3) inverse3(matrix(3) A, double& detA);