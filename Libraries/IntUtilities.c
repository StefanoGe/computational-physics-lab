#ifndef SG_INTEGER_UTILITIES_C
#define SG_INTEGER_UTILITIES_C

#include <stdlib.h>


typedef struct _arrayInt
{
	int * val;
	int length;
} ArrayInt;


ArrayInt allocArrI ( int length)
{
	ArrayInt array;
	array.length = length;
	array.val = (int *)malloc(length * sizeof(int));
	
	if (array.val == NULL)
		raiseErr("Errore nell'allocazione della memoria in funzione allocArrI\n");
		
	return array;
}

void freeArrI( ArrayInt array )
{
	free( array.val );
	array.val = NULL;
}

void printArrIPar (ArrayInt array)
{
	putchar('[');
	for(int i = 0; i < array.length; i ++)
		printf("%d, ", array.val[i]);
	printf("]\n");
}

#endif
