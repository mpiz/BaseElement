#pragma once

#include <math.h>

class CGM{
 public:
	 void init(int *gi_s, int *gj_s, double *di_s, double *gg_s, int n_s);
	 void solve(double *rp_s, double *&solution);

 private:

	 void make_LLT_decomposition();
	 void mul_matrix(double *f, double *&x);
	 void solve_L(double *f, double *&x);
	 void solve_LT(double *f, double *&x);
	 void solve_LLT(double *f, double *&x);
	 double dot_prod(double *a, double *b);

	 int n;
	 int *gi, *gj;
	 double *di, *gg, *rp, *r, *x0, *z, *p, *s;
	 double *L_di, *L_gg;
};





