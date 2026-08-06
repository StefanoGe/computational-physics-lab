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
} BarFit;

void interp_barf_init_points(BarFit *barf, const Array *points);

BarFit interp_barf_new_points(const Array *points);

void interp_barf_init_eq_points( BarFit *barf, double x1, double x2, int npoints );

BarFit interp_barf_new_eq_points(double x1, double x2, int npoints);

void interp_barf_add_func(BarFit *barf, const ParamFunc *f);

void interp_barf_add_fvalues(BarFit *barf, const ParamFunc *opt_f);

double interp_barf_get_value( const BarFit *barf, double x );

void interp_barf_free(BarFit *barf);
