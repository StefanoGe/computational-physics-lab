#ifndef ARRAY_H
#define ARRAY_H


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

void arr_init( Array *arr, int size );

void arr_free(Array *arr);

Array arr_asarr(double *data, int size);

Array arr_new(int size);

Array arr_map( const Array *x, ScalarFunc f );

void arr_free_many(Array **arrs, int n);

Array arr_linspace (double start, double end, int num_values );

#endif
