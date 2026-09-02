#ifndef SG_INTEGRATION_H
#define SG_INTEGRATION_H

#include "genutil.h"
#include "array.h"
#include "matrix.h"

#define DEF_TOL -1

typedef enum {
	BOUNDED,
	REAL_LINE,
	HALF_LINE,
	HALF_LINE_EXP
}Dbexp_mode;

double intg_trap( const ParamFuncPtr f, void *params,
			double x1, double x2, int n_subint );

double intg_simp( const ParamFuncPtr f, void *params,
					double x1, double x2, int n_subint );

double intg_legendre_root( int degree, int root_index, Array *debug );

double intg_gauss_legendre_weight( double x, int n );

double intg_legendre( double x, int degree );

double intg_legendre_wrap(double x, void *params);

double intg_gauss_legendre( ParamFuncPtr f, void *params, int order, 
			double x1, double x2 );

double intg_double_exp_real_line(ParamFuncPtr f, void *params, int N);

double intg_double_exp(ParamFuncPtr f, void *params, int N, Dbexp_mode mode,
	double a, double b,double tol);


#endif
