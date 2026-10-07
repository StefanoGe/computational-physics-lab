#ifndef ROOTF_H
#define ROOTF_H
#include "array.h"
#include "genutil.h"

#define DEF_TOL -1

double root_bis( ParamFuncPtr fnc, void *params, 
			double x1, double x2, double tol, Array *debug);

double root_newt( ParamFuncPtr f, void *params_f, 
				  ParamFuncPtr derivative, void *params_der,
				  double start, double xtol, double ftol, Array *debug );


double root_sec( const ParamFuncPtr f, void *params,
					double x1, double x2, 
					double xtol, double ftol, Array *debug, int *exit_status);

double root_newt_bracket( ParamFuncPtr f, void *params_f, 
				  ParamFuncPtr derivative, void *params_der,
				  double x1,double x2, double xtol, double ftol, Array debug[]);

double root_inv_quad_int(const ParamFuncPtr f, void *params, double x1,
				double x2, double x3, double xtol, double ftol, Array *debug,
				int *exit_status);

#endif
