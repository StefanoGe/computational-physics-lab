#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "genutil.h"
#include "array.h"
#include "matrix.h"
#include <float.h>

#define SINGULAR 0
#define NONSINGULAR 1


int interp_polyn_vand_lup(const Array *x, const Array *y, Array *coeffs,
								double *cond_number);

typedef struct {
	int npoints;
	double *points;
	double *weights;
	double *f_values;
	ParamFuncPtr func;
	void *params;
	bool are_points_projection;
	FuncPtr inv_projection;
} BarFit;

void interp_barf_init_points(BarFit *barf, const Array *points);

BarFit interp_barf_new_points(const Array *points);

void interp_barf_init_eq_points( BarFit *barf, double x1, double x2, int npoints );

BarFit interp_barf_new_eq_points(double x1, double x2, int npoints);

void interp_barf_add_func(BarFit *barf, const ParamFunc *f);

void interp_barf_add_fvalues(BarFit *barf, const ParamFunc *opt_f);

double interp_barf_get_value( const BarFit *barf, double x );

void interp_barf_free(BarFit *barf);

void interp_barf_init_cheb2_points( BarFit *barf, double x1, double x2, int npoints );

BarFit interp_barf_new_cheb2_points(double x1, double x2, int npoints);

double interp_barf_get_value_wrap(double x, void *barf);

void interp_barf_init_cheb1_points( BarFit *barf, double x1, double x2, int npoints );

BarFit interp_barf_new_cheb1_points(double x1, double x2, int npoints);

void interp_barf_add_fvalues_realline(BarFit *barf, 
	const ParamFunc *opt_f, FuncPtr projection, FuncPtr inv_projection);

double interp_real_line_to_unit( double z );

double interp_unit_to_real_line( double x );

#endif
