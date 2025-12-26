// Es 3_1.c

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "genutil.h"

/*
 * Computes the average of a sample. Float version.
 */

float aveF( float *s, int size )
{
	float sum = 0;
	for (int i = 0; i < size; i++)
		sum += s[i];
	
	return sum / (float) size;
}

float naive_varianceF( float *s, int size )
{
	if (size == 1)
		return 0;
	
	float ave = aveF( s, size );
	float sum = 0;
	for( int i = 0; i < size; i++ )
		sum += ( s[i] - ave )*( s[i] - ave );
		
	return sum / (float)(size - 1);
}

float smarter_varianceF( float *s, int size )
{
	float sum_squares = 0;
	float sum = 0;
	
	for (int i = 0; i < size; i++)
	{
		sum_squares += s[i] * s[i];
		sum += s[i];
	}
	
	eprint("sum = %f\n", sum);
	eprint("sum squares = %f\n", sum_squares);
	
	return ( sum_squares - sum*sum / (float)size ) / (float) (size - 1);
}

void testVar(float *s, int size)
{
	printf( "Bad variance is %f\n", naive_varianceF( s, size ) );
	printf( "Better Variance is %f\n", smarter_varianceF(s, size) );
}



int main()
{
	

	
	exit(EXIT_SUCCESS);
}
/*
void testAve()
{

	
	ArrayFloat testArr1 = buildArrF(3, 1,2.1f,3.4f);
	printArrFPar( testArr1 );
	printf( "Average is %f\n", aveArrF( testArr1 ) );
	printf("Hello world!");
	
	freeArrF( testArr1 );
	
}

void testVar1()
{
	ArrayFloat testArr = allocArrF(10);
	setValueArrF( testArr, 1 );

	
	ArrayFloat testArr1 = buildArrF(3, 1.0f, 2.0f, 3.0f);
	printArrFPar( testArr1 );
	printf( "Variance is %f\n", badVariance( testArr1 ) );
	
	printf( "Better Variance is %f\n", betterVariance(testArr) );
	printf( "Better Variance is %f\n", betterVariance( testArr1 ) );
	
	freeArrF( testArr );
	freeArrF( testArr1 );
}

void testSample()
{
	ArrayFloat testArr1 = buildArrF(3, 1e3, 1+1e3, 2+1e3 );
	ArrayFloat testArr2 = buildArrF(3, 1e6, 1+1e6, 2+1e6 );
	ArrayFloat testArr3 = buildArrF(3, 1e7, 1+1e7, 2+1e7 );
	ArrayFloat testArr4 = buildArrF(3, 1e8, 1+1e8, 2+1e8 );
	
	testVar(testArr1);
	testVar(testArr2);
	testVar(testArr3);
	testVar(testArr4);
	
	freeArrF(testArr1);
	freeArrF(testArr2);
	freeArrF(testArr3);
	freeArrF(testArr4);
}
*/
