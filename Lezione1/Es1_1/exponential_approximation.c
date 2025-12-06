// Prova

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comp_physics.h"

#define MAX_N 10
#define N_INTERVALS 100
#define MAX_FACTORIAL 1000

int factorial( int n )
{
	static int stored_values[MAX_FACTORIAL] = {1};
	static int count = 2;
	
	if(n >= MAX_FACTORIAL)
		raiseErr("Maximum that can be evaluated is set to %d", MAX_FACTORIAL);
	
	if ( n >= count)
		for(int i = count; i<=n; i++ )
			stored_values[i] = stored_values[i-1]*i;
	
	return stored_values[n];
}

double basicExp( double x, int N )
{
	double sum = 0;
	
	for(int i=0; i<=N; i++)
		sum += pown(x, i) / factorial(i);
		
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


/*
MatrixDouble calcAbsErrMat( int max_N, int n_intervals )
{
	const int n_points = n_intervals + 1;
	MatrixDouble absErrMat = allocMatD( n_points, max_N + 1 );
	
	ArrayDouble xAxis = linspaceD( 0, 1, n_intervals + 1);
	cpArrayToColMatD(xAxis, absErrMat, 0, true );
	
	for(int N = 1; N < max_N; N ++)
		for( int i = 0; i < n_points; i++ )
			relErrValues.val[i][j] = powInt( x, N + 1 )/fact[i +1];
	
	return absErrMat;
}
*/

// Nota: questo algoritmo è molto inefficiente. Si può fare molto di meglio.
void analysis ( int max_N, int n_intervals )
{
	const int n_xvalues = n_intervals + 1;
	MatrixDouble absErrMatrix = allocMatD( max_N + 1, n_xvalues );
	ArrayDouble xAxis = linspaceD( 0, 1, n_xvalues);
	// To the first row we assign the values of x axis for printing.
	passArrayToRowMatD(absErrMatrix, 0, xAxis);
	
	// We load the matrix with abs errs values.
	for(int i = 1; i <= max_N; i++)
		passArrayToRowMatD( absErrMatrix, i, 
							absErrArr( i, n_xvalues) );
	
	MatrixDouble relErrValues = allocMatD( max_N + 1, n_xvalues );
	cpArrayToRowMatD(relErrValues, 0, xAxis);
	
	// We compute the relative error matrix.
	for( int i = 1; i <= max_N; i++ )
		for(int j = 1; j < n_xvalues; j++)
			relErrValues.val[i][j] 
			= absErrMatrix.val[i][j]*factorial(i +1)
			/pown( absErrMatrix.val[0][j], i + 1 );
	
	printDatMatD( absErrMatrix, "e_approx.dat",   "%.14lf ", true);
	printDatMatD( relErrValues, "rel_errors.dat", "%.14lf ", true);
	
	freeAllMatD(absErrMatrix, relErrValues, NULL_MAT);
}


int main()
{
	analysis(MAX_N, N_INTERVALS);
	
	return 0;
}
