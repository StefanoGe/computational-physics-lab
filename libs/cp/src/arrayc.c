#include <stdlib.h>
#include "arrayc.h"
#include "genutil.h"
#include <string.h>

void arrc_init(ArrayC *arr, int size)
{
	if(arr->owns_data)
		SAFE_FREE(arr->data);
	
	SAFE_ALLOC(arr->data, size);
	if(!arr->data)
		raiseErr("Error in memory allocation");
		
	arr->size = size;
	arr->capacity = size;
	arr->owns_data = true;
}

void arrc_free(ArrayC *arr)
{
	if(arr->owns_data)
		SAFE_FREE(arr->data);

	arr->data = NULL;
	arr->size = 0;
	arr->capacity = 0;
	arr->owns_data = false;
}

void arrc_reserve(ArrayC *arr, int capacity)
{
	if(capacity <= arr->capacity)
		return;

	if(!arr->owns_data && arr->data != NULL)
		raiseErr("Cannot reserve to non-owning array");
        
	arr->owns_data = true;

	double complex *new_data =
		realloc(arr->data, sizeof(double complex) * capacity);

	if(!new_data)
		raiseErr("Error in memory allocation");

	arr->data = new_data;
	arr->capacity = capacity;
}

void arrc_append(ArrayC *arr, double complex x)
{
	if(!arr->owns_data && arr->data != NULL)
		raiseErr("Cannot append to non-owning array");
	
	if(arr->size >= arr->capacity)
	{
		int new_capacity = (arr->capacity == 0 ? 1 : arr->capacity * 2);
		arrc_reserve(arr, new_capacity);
	}

	arr->data[arr->size] = x;
	arr->size++;
}

ArrayC arrc_asarr(double complex *data, int size)
{
	ArrayC arr;
	arr.data = data;
	arr.size = size;
	arr.owns_data = false;
	arr.capacity = size;

	return arr;
}

ArrayC arrc_new(int size)
{
	ArrayC arr = {0};
	arrc_init(&arr, size);

	return arr;
}

ArrayC arrc_map(const ArrayC *x, FuncPtrC f)
{
	ArrayC arr = arrc_new(x->size);

	for(int i = 0; i < x->size; i++)
		ARR(arr,i) = f(ARRP(x,i));

	return arr;
}

ArrayC arrc_map_par(const ArrayC *x, ParamFuncPtrC f, void *params)
{
	ArrayC arr = arrc_new(x->size);

	for(int i = 0; i < x->size; i++)
		ARR(arr,i) = f(ARRP(x,i), params);

	return arr;
}

void arrc_free_many(ArrayC **arrs, int n)
{
	for(int i = 0; i < n; i++)
		arrc_free(arrs[i]);
}

ArrayC arrc_linspace(double complex start, double complex end, int num_values)
{
	if(num_values < 2)
		raiseErr("You must ask for at least 2 values in function arrc_linspace");
		
	ArrayC arr = arrc_new(num_values);
	
	const int n_increments = num_values - 1;
	const double complex step_size = (end - start) / n_increments;
	
	for(int i = 0; i < num_values; i++)
		ARR(arr,i) = step_size * i + start;
	
	return arr;
}

void arrc_print(FILE *file, const ArrayC *arr, const char *format)
{
	for(int i = 0; i < arr->size; i++)
		fprintf(file, format, ARRP(arr,i));
}

void arrc_print_inline(const ArrayC *arr, const char *format, bool parentheses,
	bool linebreak)
{
	if(parentheses)
		fputs("[ ", stdout);
	
	for(int i = 0; i < arr->size; i++)
	{
		if(i != 0)
			fputs(", ", stdout);

		printf(format, ARRP(arr, i));
	}
	
	if(parentheses)
		fputs(" ]", stdout);

	if(linebreak)
		putchar('\n');
}

void arrc_cp(const ArrayC *src, ArrayC *dest)
{
	arrc_init(dest, src->size);
	memcpy(dest->data, src->data, src->size * sizeof *dest->data);
}

void arrc_setv(ArrayC *arr, double complex v)
{
	for(int i = 0; i < arr->size; i++)
		ARRP(arr,i) = v;
}

void arrc_scale(ArrayC *x, double complex alpha)
{
	for(int i = 0; i < x->size; i++)
		ARRP(x,i) *= alpha;
}

void arrc_axpy(ArrayC *y, double complex alpha, const ArrayC *x)
{
	if(x->size != y->size)
		raiseErr("arrays should be of same dimension");

	for(int i = 0; i < x->size; i++)
		ARRP(y,i) += alpha * ARRP(x,i);
}

void arrc_null(ArrayC *arr)
{
	arr->data = NULL;
	arr->size = 0;
	arr->capacity = 0;
	arr->owns_data = false;
}

ArrayC arrc_slice(const ArrayC *arr, int start, int end)
{
	if(start >= end)
		raiseErr("end must be strictly higher than start: "
				"start: %d end: %d", start, end);
	
	if(end > arr->size)
		raiseErr("end must be less or equal to arr->size: "
				"arr->size: %d end: %d", arr->size, end);
	
	const int new_size = end - start;
	ArrayC new = arrc_new(new_size);

	for(int i = 0; i < new_size; i++)
		ARR(new,i) = ARRP(arr,i + start);

	return new;
}

Array arrc_real(const ArrayC *arr)
{
	Array new=arr_new(arr->size);
	for(int i=0; i<arr->size; i++)
		ARR(new,i)=creal(ARRP(arr,i));
	return new;
}
