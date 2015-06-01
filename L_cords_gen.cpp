#include "L_cords_gen.h"


double determenant3(matrix(3) A) {

	return A[0][0]*A[1][1]*A[2][2] + A[1][0]*A[2][1]*A[0][2] + A[0][1]*A[1][2]*A[2][0] - A[2][0]*A[1][1]*A[0][2] - A[2][1]*A[1][2]*A[0][0] - A[1][0]*A[0][1]*A[2][2];
}

matrix(4) inverse4(matrix(4) A, double& detA) {
	matrix(3) B; //вспомогатльная матрица
	matrix(4) Inv_A; //обращенная матрица
	
	detA = 0;
	for(int a_j = 0; a_j < 4; a_j++){

		for(int i = 1; i < 4; i++)
			for(int j = 0; j < a_j; j++)
				B[i-1][j] = A[i][j];

		for(int i = 1; i < 4; i++)
			for(int j = a_j+1; j < 4; j++)
				B[i-1][j-1] = A[i][j];
		double sum_el = A[0][a_j]*determenant3(B);
		if((a_j + 2)%2 == 1)
			sum_el *= -1;
		detA += sum_el;
	}
	
	for(int inv_j = 0; inv_j < 4; inv_j++) {
		for(int inv_i = 0; inv_i < 4; inv_i++) {
			//Слева, сверху от элемента
			for(int i = 0; i < inv_i; i++)
				for(int j = 0; j < inv_j; j++)
					B[i][j] = A[i][j];
			//Слева, снизу от элемента
			for(int i = inv_i+1; i < 4; i++)
				for(int j = 0; j < inv_j; j++)
					B[i-1][j] = A[i][j];
			//Справа, сверху
			for(int i = 0; i < inv_i; i++)
				for(int j = inv_j+1; j < 4; j++)
					B[i][j-1] = A[i][j];
			//Справа, снизу
			for(int i = inv_i+1; i < 4; i++)
				for(int j = inv_j+1; j < 4; j++)
					B[i-1][j-1] = A[i][j];
			
			Inv_A[inv_j][inv_i] = determenant3(B) / detA;
			if((inv_i + inv_j + 2)%2 == 1)
				Inv_A[inv_j][inv_i] *= -1;
		}
	}
	return Inv_A;
}

matrix(3) inverse3(matrix(3) A, double& detA) {
	matrix(3) Inv_A;
	detA = determenant3(A);

	Inv_A[0][0] = A[1][1]*A[2][2] - A[2][1]*A[1][2];
	Inv_A[0][1] = -(A[0][1]*A[2][2] - A[2][1]*A[0][2]);
	Inv_A[0][2] = A[0][1]*A[1][2] - A[1][1]*A[0][2];

	Inv_A[1][0] = -(A[1][0]*A[2][2] - A[2][0]*A[1][2]);
	Inv_A[1][1] = A[0][0]*A[2][2] - A[2][0]*A[0][2];
	Inv_A[1][2] = -(A[0][0]*A[1][2] - A[1][0]*A[0][2]);

	Inv_A[2][0] = A[1][0]*A[2][1] - A[2][0]*A[1][1];
	Inv_A[2][1] = -(A[0][0]*A[2][1] - A[2][0]*A[0][1]);
	Inv_A[2][2] = A[0][0]*A[1][1] - A[1][0]*A[0][1];

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			Inv_A[i][j] /= detA;

	return Inv_A;
}