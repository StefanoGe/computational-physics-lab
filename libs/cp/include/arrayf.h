#ifndef ARRAYF_H
#define ARRAYF_H

#include <stdio.h>
#include <stdbool.h>
#include "genutil.h"

// -----------------
#ifndef ARR_ACCESS
	#define ARR_ACCESS

	#define ARR(A,i) ((A).data[(i)])
	#define ARRP(A,i) ((A)->data[(i)])
#endif
// -----------------

typedef float(FuncPtrF)(float);

typedef float(ParamFuncPtrF)(float,void*);

/*
data == NULL && size == 0 && capacity == 0
    => empty, appendable ArrayF

data != NULL && owns_data == false
    => non-owning view, cannot grow

data != NULL && owns_data == true
    => owning, growable ArrayF
*/
typedef struct {
	float *data;
	int size;
	int capacity;
	bool owns_data;
} ArrayF;

void farr_init(ArrayF *, int size);

void farr_free(ArrayF *);

ArrayF farr_asarr(float *data, int size);

ArrayF farr_new(int size);

ArrayF farr_map(const ArrayF *, FuncPtrF f);

ArrayF farr_map_par(const ArrayF *x, ParamFuncPtrF f, void *params);

void farr_free_many(ArrayF **arrs, int n);

ArrayF farr_linspace(float start, float end, int num_values);

void farr_print(FILE *, const ArrayF *, const char *format);

void farr_print_inline(const ArrayF *arr, const char *format, bool parentheses,
		bool linebreak);

void farr_cp(const ArrayF *src, ArrayF *dest);

void farr_setv(ArrayF *arr, float v);

void farr_axpy(ArrayF *y, float alpha, const ArrayF *x);

void farr_scale(ArrayF *x, float alpha);

void farr_append(ArrayF *arr, float x);

void farr_reserve(ArrayF *arr, int capacity);

void farr_null(ArrayF *arr);

ArrayF farr_slice(const ArrayF *arr, int start, int end);

float farr_max(const ArrayF *arr);

#endif
