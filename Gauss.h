#pragma once

#include <math.h>
#include "macro.h"

template<int n> bool SLAE_solution_Gauss(matrix(n)& A, array<double, n>& b, array<double, n>& x); //решение СЛАУ, метод Гаусса-Жордана

template<int n> bool trianglematrix1(matrix(n)& A, array<double, n>& x); //Приведение матрицы к верхнему треугольному виду
template<int n> void transf1(matrix(n)& A, array<double, n>& x, int i); // перестановка строк


// ====================== Реализации =========================

#include "Gauss.h"

template<int n> bool SLAE_solution_Gauss(matrix(n)& A, array<double, n>& b, array<double, n>& x){
	bool bad;
	for (int i = 0; i < n; i++)
		x[i] = b[i];
	bad = trianglematrix1(A, x);
	if (!bad) return bad;
	for (int i = n - 1; i >= 0; i--){
		double sum = 0;
		for (int j = i + 1; j < n; j++)
			sum += A[i][j] * x[j];
		x[i] -= sum;
	}
	return true;
}


template<int n> bool trianglematrix1(matrix(n)& A, array<double, n>& x){
	for (int i = 0; i < n; i++){
		transf1(A, x, i);
		double A_d = A[i][i];
		if (fabs(A_d) < 1E-20) return false;

		for (int p = i; p < n; p++){
			A[i][p] /= A_d;
		}
		x[i] /= A_d;

		for (int j = i + 1; j < n; j++){
			double A_j = A[j][i];
			if (A_j != 0){
				for (int k = i; k < n; k++)
					A[j][k] -= A[i][k] * A_j;
				x[j] -= x[i] * A_j;

			}
		}
	}
	return true;
}


template<int n> void transf1(matrix(n)& A, array<double, n>& x, int i){

	int line = i;
	for (int j = i + 1; j < n; j++)
		if (fabs(A[j][i]) > fabs(A[line][i]))
			line = j;

	if (line != i){
		double c;
		for (int j = i; j < n; j++){
			c = A[i][j];
			A[i][j] = A[line][j];
			A[line][j] = c;
		}
		c = x[i];
		x[i] = x[line];
		x[line] = c;
	}
}


