#ifndef ARRAYF_H
#define ARRAYF_H

#include <stdio.h>
#include "array.h"

#ifndef ARR_ACCESS
#define ARR_ACCESS

	#define ARR(A,i) ((A).data[(i)])
	#define ARRP(A,i) ((A)->data[(i)])

#endif

typedef struct {
	float *data;
	int size;
	bool owns_data;
} ArrayF;

void arrf_init( ArrayF *arr, int size );

void arrf_free(ArrayF *arr);

ArrayF arrf_asarr(float *data, int size);

ArrayF arrf_new(int size);

void arrf_print( FILE *file, const ArrayF *arr, const char *format );

Array arr_from_float( ArrayF *arrf);

// void arrf_free_many(ArrayF **arrs, int n);

// Array arr_linspace (double start, double end, int num_values );

#endif
