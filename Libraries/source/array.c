#include <stdlib.h>
#include "array.h"
#include "genutil.h"

void arr_init( Array *arr, int size )
{
	if(arr->owns_data)
		free(arr->data);
	
	arr->data = malloc(sizeof(double) * size);
	if( !arr->data )
		raiseErr("Error in memory allocation");
		
	arr->size = size;
	arr->owns_data=true;
}

void arr_free(Array *arr)
{
	if(arr->owns_data)
		free( arr->data );
	arr->data=nullptr;
	arr->size=0;
	arr->owns_data=false;
}

Array arr_asarr(double *data, int size)
{
	Array arr;
	arr.data=data;
	arr.size=size;
	return arr;
}

Array arr_new(int size)
{
	Array arr={0};
	arr_init(&arr, size);
	return arr;
}

Array arr_map( const Array *x, ScalarFunc f )
{
	Array arr=arr_new(x->size);
	for(int i=0; i<x->size; i++)
		ARR(arr,i)=f(ARRP(x,i));
	return arr;
}

void arr_free_many(Array **arrs, int n)
{
    for (int i = 0; i < n; i++)
        arr_free(arrs[i]);
}

Array arr_linspace (double start, double end, int num_values )
{
	if (num_values < 2)
		raiseErr("You must ask for at least 2 values in function linspaceD");
		
	Array arr = arr_new(num_values);
	
	const int n_increments = num_values - 1;
	const double step_size = (end - start)/n_increments;
	
	for(int i=0; i < num_values; i++)
		ARR(arr,i) = step_size*i + start;
	
	return arr;
}

