// Float.c

#ifndef SG_FLOAT_UTILITIES_C
#define SG_FLOAT_UTILITIES_C

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

typedef struct _arrayFloat
{
	float * val;
	int length;
	
} ArrayFloat;

ArrayFloat allocArrF ( int length)
{
	ArrayFloat array;
	array.length = length;
	array.val = (float *)malloc(length * sizeof(float));
	if (array.val == NULL) {
		printf("Errore nell'allocazione della memoria\n");
		exit(1);
    }
	return array;
}

ArrayFloat buildArrF ( int length, ... )
{
	va_list list;
	va_start(list, length);
	
	ArrayFloat array = allocArrF( length );
	
	for( int i = 0; i < length; i++ )
		array.val[ i ] = (float)va_arg( list, double );
		
	va_end(list);
	
	return array;
}

void freeArrF( ArrayFloat array )
{
	free( array.val );
	array.val = NULL;
}

void printArrFPar (ArrayFloat array)
{
	putchar('[');
	for(int i = 0; i < array.length; i ++)
		printf("%.10f, ", array.val[i]);
	printf("]\n");
}

void setValueArrF(ArrayFloat array, float value)
{
	for( int i = 0; i < array.length; i++ )
		array.val[i] = value;
}


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

#endif

