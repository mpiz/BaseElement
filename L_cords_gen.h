#pragma once
//Построение L-координат для треугольника и тетраэдра

#include "macro.h"

//========= Тетраэдры ==========

double determenant3(matrix(3) A); //определетель матрицы 3x3
matrix(4) inverse4(matrix(4) A, double& detA); //обращение матрицы 4x4 (метод Крамера)

//========= Треугольники ==========

matrix(3) inverse3(matrix(3) A, double& detA);