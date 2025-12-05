// Es 3_1double.c

#include "comp_physics.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


double aveArrD( ArrayDouble array )
{
	double sum = 0;
	for (int i = 0; i < array.length; i++)
		sum += array.val[i];
	
	return sum / (double) array.length;
}

double badVariance( ArrayDouble array )
{
	if (array.length == 1)
		return 0;
	
	double ave = aveArrD( array );
	double sum = 0;
	for( int i = 0; i < array.length; i++ )
		sum += ( array.val[i] - ave )*( array.val[i] - ave );
		
	return sum / (double) (array.length - 1);
}

double betterVariance( ArrayDouble array )
{
	double sum_squares = 0;
	double sum = 0;
	
	for (int i = 0; i < array.length; i++)
	{
		sum_squares += array.val[i] * array.val[i];
		sum += array.val[i];
	}
	
	printf("sum = %f\n", sum);
	printf("sum squares = %f\n", sum_squares);
	
	return ( sum_squares - sum*sum / (double)array.length ) / 
		(double) (array.length - 1);
}

void testVar(ArrayDouble array)
{
	printf( "Bad variance is %lf\n", badVariance( array ) );
	printf( "Better Variance is %lf\n", betterVariance(array) );
}

void testSample()
{
	ArrayDouble testArr1 = buildArrD(3, 1e3, 1+1e3, 2+1e3 );
	ArrayDouble testArr2 = buildArrD(3, 1e6, 1+1e6, 2+1e6 );
	ArrayDouble testArr3 = buildArrD(3, 1e7, 1+1e7, 2+1e7 );
	ArrayDouble testArr4 = buildArrD(3, 1e8, 1+1e8, 2+1e8 );
	
	testVar(testArr1);
	testVar(testArr2);
	testVar(testArr3);
	testVar(testArr4);
	
	freeArrD(testArr1);
	freeArrD(testArr2);
	freeArrD(testArr3);
	freeArrD(testArr4);
}

int main()
{
	
	testSample();
	
	exit(EXIT_SUCCESS);
}
/*
void testAve()
{

	
	ArrayDouble testArr1 = buildArrD(3, 1.0,2.1,3.4);
	printArrDPar( testArr1 );
	printf( "Average is %lf\n", aveArrD( testArr1 ) );
	printf("Hello world!");
	
	freeArrD( testArr1 );
	
}

void testVar1()
{
	ArrayDouble testArr = allocArrD(10);
	setValueArrD( testArr, 1 );

	
	ArrayDouble testArr1 = buildArrD(3, 1.0, 2.0, 3.0);
	printArrDPar( testArr1 );
	printf( "Variance is %lf\n", badVariance( testArr1 ) );
	
	printf( "Better Variance is %lf\n", betterVariance(testArr) );
	printf( "Better Variance is %lf\n", betterVariance( testArr1 ) );
	
	freeArrD( testArr );
	freeArrD( testArr1 );
}

*/
