// Prova

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "DoubleUtilities.c"

int fact( int n )
{
	if (n==0)
		return 1;
		
	return fact(n-1) * n;
}

double powInt( double x, int n)
{
	if (n == 0)
		return 1;
		
	return x * powInt(x, n -1);
}

double basicExp( double x, int N )
{
	double sum = 0;
	
	for(int i=0; i<=N; i++)
		sum += powInt(x, i) / fact(i);
		
	return sum;
}

double dankExp( double x, int N )
{
	double sum = 1;
	for( int i = 0; i < N; i++ )
	{
		sum *= x;
		sum /= ( N - i );
		sum += 1;
	}
	return sum;
}

double absoluteError( double x, int N )
{
	return fabs( dankExp( x, N ) - exp( x ) );
}

ArrayDouble absErrArr( int n, int length )
{
	ArrayDouble array = allocArrD( length);
	for (int i = 0; i < length; i ++)
		array.val[i] = absoluteError( (double) i / ( length - 1), n );
	return array;
}
// Nota: questo algoritmo è molto inefficiente. Si può fare molto di meglio.
void analysis ( int max_N, int n_intervals )
{
	MatrixDouble matrix = allocMatD( max_N + 1, n_intervals + 1 );
	ArrayDouble xAxis = linspaceD( 0, 1, n_intervals + 1);
	passArrayToRowMatD(matrix, 0, xAxis);
	
	for(int i = 1; i <= max_N; i++)
		passArrayToRowMatD( matrix, i, 
							absErrArr( i, n_intervals + 1) );
	
	
	MatrixDouble traMatrix = transposeMatD(matrix);
	
	printDatMatD(traMatrix, "e_approx.dat", "%.14lf ");
	
	freeMatD(matrix);
	freeMatD(traMatrix);
}


int main()
{
	printf("Hello world!");
	analysis(10, 100);
	
	return 0;
}



/*
 *	ArrayFloat arr = allocArrF(10);
 *	setValueArrF(arr, 5);
 *	printArrFPar(arr);
 *	printf("Hello world!\n");
 */

//void testExps()
//{
	//printf( "%lf\n", basicExp(1, 5) );
	//printf( "%lf\n", dankExp(1, 5) );
//}
