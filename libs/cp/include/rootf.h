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



#endif
