#include "genutil.h"
#include "array.h"
#include "matrix.h"

#define SINGULAR 0
#define NONSINGULAR 1


int interp_polyn_vand_lup(const Array *x, const Array *y, Array *coeffs,
								double *cond_number);
