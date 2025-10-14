// Esercizio 1 della lezione 2

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846

typedef struct _arrayFloat
{
	float * val;
	int length;
	
} arrayFloat;

arrayFloat allocArrF ( int length)
{
	arrayFloat array;
	array.length = length;
	array.val = (float *)malloc(length * sizeof(float));
	return array;
}

void freeArrF( arrayFloat array )
{
	free( array.val );
}

float sumSquaresFloat( int maxIndex )
{
	float sum = 0;
	for ( int i = 1; i <= maxIndex; i ++ )
		sum += 1.f / (float) (i * i);
		
	return sum;
}

float revSumSquaresFloat( int maxIndex )
{
	float sum = 0;
	for ( int i = maxIndex; i > 0; i -- )
		sum += 1.f / (float) (i * i);
		
	return sum;
}

/*
void sumSquaresAnalysis( int maxIndex )
{
	
}

*/



int main ()
{


	
}

void testSumSquare()
{
	printf("Ordine normale:\n");
	
	printf("%f\n", sumSquaresFloat( 30 ) );
	
	printf("Ordine reversed:\n");
	
	printf("%f\n", revSumSquaresFloat( 30 ) );
	
	printf( "Valore vero:\n%f", PI*PI / 6 );
}
