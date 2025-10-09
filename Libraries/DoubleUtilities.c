// Double.c

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

typedef struct _arrayDouble
{
	double * val;
	int length;
	
} ArrayDouble;

ArrayDouble allocArrD ( int length)
{
	ArrayDouble array;
	array.length = length;
	array.val = (double *)malloc(length * sizeof(double));
	if (array.val == NULL) {
		printf("Errore nell'allocazione della memoria\n");
		exit(1);
    }
	return array;
}

ArrayDouble buildArrD ( int length, ... )
{
	va_list list;
	va_start(list, length);
	
	ArrayDouble array = allocArrD( length );
	
	for( int i = 0; i < length; i++ )
		array.val[ i ] = va_arg( list, double );
		
	va_end(list);
	
	return array;
}

void freeArrD( ArrayDouble array )
{
	free( array.val );
	array.val = NULL;
}

void printArrDPar (ArrayDouble array)
{
	putchar('[');
	for(int i = 0; i < array.length; i ++)
		printf("%.10lf, ", array.val[i]);
	printf("]\n");
}

void setValueArrD(ArrayDouble array, double value)
{
	for( int i = 0; i < array.length; i++ )
		array.val[i] = value;
}


ArrayDouble diffArrD( ArrayDouble array1, ArrayDouble array2)
{
	if (array1.length != array2.length )
	{
		printf("Tried to subtract arrays of different sizes in function "
		"diffArrD\n");
		exit(EXIT_FAILURE);
	}
	
	ArrayDouble diffArray = allocArrD( array1.length );
	
	for ( int i = 0; i < array1.length; i++ )
		diffArray.val[i] = array1.val[i] - array2.val[i];
		
	return diffArray;
}



