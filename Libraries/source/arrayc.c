#include "arrayc.h"
#include <stdlib.h>
#include "comp_physics.h"
#include <string.h>

void arrc_init( ArrayC *arr, int size )
{
	if(arr->owns_data)
		free(arr->data);
	
	arr->data = malloc(sizeof(double complex) * size);
	if( !arr->data )
		raiseErr("Error in memory allocation");
		
	arr->size = size;
	arr->owns_data=true;
}

void arrc_free(ArrayC *arr)
{
	if(arr->owns_data)
		free( arr->data );
	arr->data=nullptr;
	arr->size=0;
	arr->owns_data=false;
}

void arrc_free_many(ArrayC **arrs, int n)
{
    for (int i = 0; i < n; i++)
        arrc_free(arrs[i]);
}

ArrayC arrc_new(int size)
{
	ArrayC arr={0};
	arrc_init(&arr, size);
	return arr;
}

void arrc_scale(ArrayC *dest, double complex x)
{
	for(int i=0; i<dest->size; i++)
		ARRP(dest, i)*=x;
}

void arrc_add_scaled(ArrayC *dest, double complex x, const ArrayC *src)
{
    if (dest->size != src->size)
        raiseErr("Array size mismatch");

    for (int i = 0; i < dest->size; i++)
        ARRP(dest,i) += x * ARRP(src,i);
}

void arrc_cp(double complex *src, int size, ArrayC *dest)
{
	if( !dest->owns_data || dest->size!=size )
		arrc_init(dest, size);
	memcpy(dest->data, src, size*sizeof(double complex));
}





