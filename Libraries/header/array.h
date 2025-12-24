#ifndef ARRAY_H
#define ARRAY_H

#include <stdio.h>

#ifndef ARR_ACCESS
#define ARR_ACCESS

#define ARR(A,i) ((A).data[(i)])
#define ARRP(A,i) ((A)->data[(i)])

#endif

typedef double(ScalarFunc)(double);

typedef struct {
	double *data;
	int size;
	bool owns_data;
} Array;

void arr_init( Array *, int size );

void arr_free(Array *);

Array arr_asarr(double *data, int size);

Array arr_new(int size);

Array arr_map( const Array *, ScalarFunc f );

void arr_free_many(Array **arrs, int n);

Array arr_linspace (double start, double end, int num_values );

void arr_print( FILE *, const Array *, const char *format );

#endif
