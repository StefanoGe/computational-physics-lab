#ifndef ARRAYC_H
#define ARRAYC_H

#include <complex.h>

#ifndef ARR_ACCESS
#define ARR_ACCESS

#define ARR(A,i) ((A).data[(i)])
#define ARRP(A,i) ((A)->data[(i)])

#endif

typedef struct {
	double complex *data;
	int size;
	bool owns_data;
} ArrayC;

void arrc_init( ArrayC *arr, int size );

void arrc_free(ArrayC *arr);

void arrc_free_many(ArrayC **arrs, int n);

ArrayC arrc_new(int size);

void arrc_cp(double complex *src, int size, ArrayC *dest);

void arrc_add_scaled(ArrayC *dest, double complex x, const ArrayC *src);

void arrc_scale(ArrayC *dest, double complex x);

#endif
