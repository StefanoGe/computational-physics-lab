// Es 3_1.c

#include "comp_physics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


float aveArrF( ArrayFloat array )
{
	float sum = 0;
	for (int i = 0; i < array.length; i++)
		sum += array.val[i];
	
	return sum / (float) array.length;
}

float badVariance( ArrayFloat array )
{
	if (array.length == 1)
		return 0;
	
	float ave = aveArrF( array );
	float sum = 0;
	for( int i = 0; i < array.length; i++ )
		sum += ( array.val[i] - ave )*( array.val[i] - ave );
		
	return sum / (float) (array.length - 1);
}

float betterVariance( ArrayFloat array )
{
	float sum_squares = 0;
	float sum = 0;
	
	for (int i = 0; i < array.length; i++)
	{
		sum_squares += array.val[i] * array.val[i];
		sum += array.val[i];
	}
	
	printf("sum = %f\n", sum);
	printf("sum squares = %f\n", sum_squares);
	
	return ( sum_squares - sum*sum / (float)array.length ) / 
		(float) (array.length - 1);
}

void testVar(ArrayFloat array)
{
	printf( "Bad variance is %f\n", badVariance( array ) );
	printf( "Better Variance is %f\n", betterVariance(array) );
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

int main()
{
	
	testSample();
	
	exit(EXIT_SUCCESS);
}

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

