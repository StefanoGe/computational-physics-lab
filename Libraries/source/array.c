#include <stdlib.h>
#include "array.h"
#include "comp_physics.h"

void arr_init( Array *arr, int size )
{
	if( arr->data )
		free(arr->data);
	
	arr->data = malloc(sizeof(double) * size);
	if( !arr->data )
		raiseErr("Error in memory allocation");
		
	arr->size = size;
}

void arr_free(Array *arr)
{
	free( arr->data );
	arr->data=nullptr;
	arr->size=0;
}

Array arr_asarr(double *data, int size)
{
	Array arr;
	arr.data=data;
	arr.size=size;
	return arr;
}


