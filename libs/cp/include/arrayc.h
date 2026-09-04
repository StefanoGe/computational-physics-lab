#ifndef ARRAYC_H
#define ARRAYC_H

#include <stdio.h>
#include <stdbool.h>
#include <complex.h>
#include "genutil.h"

// -----------------
#ifndef ARR_ACCESS
	#define ARR_ACCESS

	#define ARR(A,i) ((A).data[(i)])
	#define ARRP(A,i) ((A)->data[(i)])
#endif
// -----------------

typedef double complex(FuncPtrC)(double complex);
typedef double complex(ParamFuncPtrC)(double complex,void*);

/*
data == NULL && size == 0 && capacity == 0
    => empty, appendable ArrayC

data != NULL && owns_data == false
    => non-owning view, cannot grow

data != NULL && owns_data == true
    => owning, growable ArrayC
*/
typedef struct {
	double complex *data;
	int size;
	int capacity;
	bool owns_data;
} ArrayC;

void arrc_init(ArrayC *, int size);

void arrc_free(ArrayC *);

ArrayC arrc_asarr(double complex *data, int size);

ArrayC arrc_new(int size);

ArrayC arrc_map(const ArrayC *, FuncPtrC f);

ArrayC arrc_map_par(const ArrayC *x, ParamFuncPtrC f, void *params);

void arrc_free_many(ArrayC **arrs, int n);

ArrayC arrc_linspace(double complex start, double complex end, int num_values);

void arrc_print(FILE *, const ArrayC *, const char *format);

void arrc_print_inline(const ArrayC *arr, const char *format, bool parentheses,
		bool linebreak);

void arrc_cp(const ArrayC *src, ArrayC *dest);

void arrc_setv(ArrayC *arr, double complex v);

void arrc_axpy(ArrayC *y, double complex alpha, const ArrayC *x);

void arrc_scale(ArrayC *x, double complex alpha);

void arrc_append(ArrayC *arr, double complex x);

void arrc_reserve(ArrayC *arr, int capacity);

void arrc_null(ArrayC *arr);

ArrayC arrc_slice(const ArrayC *arr, int start, int end);

#endif
