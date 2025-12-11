#include <stdlib.h>
#include "array.h"
#include "comp_physics.h"

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


