// Esercizio 1 della lezione 2

#include <stdio.h>
#include <stdlib.h>

#define n 10

void printArrayInt( int * array, int size )
{
	putchar( '[' );
	for ( int i =0; i < size; i++ )
		printf( "%d, ", array[i] );
	printf( "]\n" );
	
}

int main ()
{
	printf("Hello world!");
	
	int array[n] = {0};
	
	printArrayInt( array, n );
	
}
