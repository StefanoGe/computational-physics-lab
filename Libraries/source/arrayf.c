// Float.c

#include "genutil.h"
#include "arrayf.h"
#include <stdlib.h>

void arrf_init( ArrayF *arr, int size )
{
	if(arr->owns_data)
		free(arr->data);
	
	arr->data = malloc(sizeof(float) * size);
	if( !arr->data )
		raiseErr("Error in memory allocation");
		
	arr->size = size;
	arr->owns_data=true;
}

void arrf_free(ArrayF *arr)
{
	if(arr->owns_data)
		free( arr->data );
	arr->data=nullptr;
	arr->size=0;
	arr->owns_data=false;
}

ArrayF arrf_new(int size)
{
	ArrayF arr={0};
	arrf_init(&arr, size);
	return arr;
}

ArrayF arrf_asarr(float *data, int size)
{
	ArrayF arr;
	arr.data=data;
	arr.size=size;
	arr.owns_data=false;
	return arr;
}

void arrf_print( FILE *file, const ArrayF *arr, const char *format )
{
	for(int i=0; i<arr->size; i++)
		fprintf(file, format, ARRP(arr,i));
}

void arrf_set_value(ArrayF *arr, float value)
{
	for( int i=0; i < arr->size; i++ )
		ARRP(arr,i)=value;
}

Array arr_from_float( ArrayF *arrf)
{
	Array arr=arr_new(arrf->size);
	for(int i=0; i<arrf->size; i++)
		ARR(arr,i)=ARRP(arrf,i);
	return arr;
}

/*
ArrayFloat diffArrF( ArrayFloat array1, ArrayFloat array2)
{
	if (array1.length != array2.length )
	{
		printf("Tried to subtract arrays of different sizes in function "
		"diffArrF\n");
		exit(EXIT_FAILURE);
	}
	
	ArrayFloat diffArray = allocArrF( array1.length );
	
	for ( int i = 0; i < array1.length; i++ )
		diffArray.val[i] = array1.val[i] - array2.val[i];
		
	return diffArray;
}
*/

