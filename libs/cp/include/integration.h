#ifndef SG_INTEGRATION_H
#define SG_INTEGRATION_H

#include "genutil.h"
#include "array.h"
#include "matrix.h"

double intg_trap( const ParamFuncPtr f, void *params,
			double x1, double x2, int n_subint );

double intg_simp( const ParamFuncPtr f, void *params,
					double x1, double x2, int n_subint );



#endif
