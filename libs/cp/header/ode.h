#ifndef ODE_H
#define ODE_H

#include "array.h"
#include "matrix.h"

typedef double(*Derivative) (double, double *, void *);

typedef struct{
	Derivative f;
	void *params;
	int nparams;
}ParDer;

void ode_feuler( double init_time, double final_time, const double init_value[], 
				const Derivative ders[], void **f_parameters,
				int sys_size, int nsteps, Array *t_output, Matrix *f_output );

void ode_ie2( double init_time, double final_time, const double init_value[], 
				const Derivative ders[], void **f_parameters,
				int sys_size, int nsteps, Array *t_output, Matrix *f_output );

void ode_rk4( double init_time, double final_time, const double init_value[], 
				const Derivative ders[], void **f_parameters,
				int sys_size, int nsteps, Array *t_output, Matrix *f_output );

#endif
