#include <stdlib.h>
#include "arrayf.h"
#include "genutil.h"
#include <string.h>

void farr_init( ArrayF *arr, int size )
{
	if(arr->owns_data)
		SAFE_FREE(arr->data);
	
	SAFE_ALLOC(arr->data,size);
	if( !arr->data )
		raiseErr("Error in memory allocation");
		
	arr->size = size;
	arr->capacity = size;
	arr->owns_data = true;
}

void farr_free(ArrayF *arr)
{
	if(arr->owns_data)
		SAFE_FREE( arr->data );

	arr->data = NULL;
	arr->size = 0;
	arr->capacity = 0;
	arr->owns_data = false;
}

void farr_reserve(ArrayF *arr, int capacity)
{
    if (capacity <= arr->capacity)
        return;

    if (!arr->owns_data && arr->data != NULL)
        raiseErr("Cannot reserve to non-owning array");
        
    arr->owns_data = true;

    float *new_data = realloc(arr->data, sizeof(float) * capacity);

    if (!new_data)
        raiseErr("Error in memory allocation");

    arr->data = new_data;
    arr->capacity = capacity;
}

void farr_append(ArrayF *arr, float x)
{
    if (!arr->owns_data && arr->data != NULL)
        raiseErr("Cannot append to non-owning array");
	
    if (arr->size >= arr->capacity)
    {
        int new_capacity = (arr->capacity == 0 ? 1 : arr->capacity * 2);
        farr_reserve(arr, new_capacity);
    }

    arr->data[arr->size] = x;
    arr->size++;
}

ArrayF farr_asarr(float *data, int size)
{
	ArrayF arr;
	arr.data = data;
	arr.size = size;
	arr.capacity = size;
	arr.owns_data = false;

	return arr;
}

ArrayF farr_new(int size)
{
	ArrayF arr = {0};
	farr_init(&arr, size);

	return arr;
}

ArrayF farr_map(const ArrayF *x, FuncPtrF f)
{
	ArrayF arr = farr_new(x->size);

	for(int i = 0; i < x->size; i++)
		ARR(arr,i) = f(ARRP(x,i));

	return arr;
}

ArrayF farr_map_par(const ArrayF *x, ParamFuncPtrF f, void *params)
{
	ArrayF arr = farr_new(x->size);

	for(int i = 0; i < x->size; i++)
		ARR(arr,i) = f(ARRP(x,i), params);

	return arr;
}

void farr_free_many(ArrayF **arrs, int n)
{
	for(int i = 0; i < n; i++)
		farr_free(arrs[i]);
}

ArrayF farr_linspace(float start, float end, int num_values)
{
	if(num_values < 2)
		raiseErr("You must ask for at least 2 values in function farr_linspace");

	ArrayF arr = farr_new(num_values);

	const int n_increments = num_values - 1;
	const float step_size = (end - start) / (float)n_increments;

	for(int i = 0; i < num_values; i++)
		ARR(arr,i) = step_size*(float)i + start;

	return arr;
}

void farr_print(FILE *file, const ArrayF *arr, const char *format)
{
	for(int i = 0; i < arr->size; i++)
		fprintf(file, format, ARRP(arr,i));
}

void farr_print_inline(const ArrayF *arr, const char *format, bool parentheses,
	bool linebreak)
{
	if(parentheses)
		fputs("[ ", stdout);

	for(int i = 0; i < arr->size; i++)
	{
		if(i != 0)
			fputs(", ", stdout);

		printf(format, ARRP(arr,i));
	}

	if(parentheses)
		fputs(" ]", stdout);

	if(linebreak)
		putchar('\n');
}

void farr_cp(const ArrayF *src, ArrayF *dest)
{
	farr_init(dest, src->size);
	memcpy(dest->data, src->data, src->size * sizeof *dest->data);
}

void farr_setv(ArrayF *arr, float v)
{
	for(int i = 0; i < arr->size; i++)
		ARRP(arr,i) = v;
}

void farr_scale(ArrayF *x, float alpha)
{
	for(int i = 0; i < x->size; i++)
		ARRP(x,i) *= alpha;
}

void farr_axpy(ArrayF *y, float alpha, const ArrayF *x)
{
	if(x->size != y->size)
		raiseErr("arrays should be of same dimension");

	for(int i = 0; i < x->size; i++)
		ARRP(y,i) += alpha * ARRP(x,i);
}

void farr_null(ArrayF *arr)
{
	arr->data = NULL;
	arr->size = 0;
	arr->capacity = 0;
	arr->owns_data = false;
}

ArrayF farr_slice(const ArrayF *arr, int start, int end)
{
	if(start >= end)
		raiseErr("end must be strictly higher than start: "
				"start: %d end: %d", start, end);

	if(end > arr->size)
		raiseErr("end must be less or equal to arr->size: "
				"arr->size: %d end: %d", arr->size, end);

	const int new_size = end - start;
	ArrayF new = farr_new(new_size);

	for(int i = 0; i < new_size; i++)
		ARR(new,i) = ARRP(arr,i + start);

	return new;
}

float farr_max(const ArrayF *arr)
{
	float max = ARRP(arr,0);

	for(int i = 1; i < arr->size; i++)
		if(max < ARRP(arr,i))
			max = ARRP(arr,i);

	return max;
}
