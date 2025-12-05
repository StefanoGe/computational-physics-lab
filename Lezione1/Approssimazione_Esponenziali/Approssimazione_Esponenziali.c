// Prova

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comp_physics.h"

#define MAX_N 10
#define N_INTERVALS 100

int fact[MAX_N * 2];

int factorial( int n )
{
	if (n==0)
		return 1;
		
	return factorial(n-1) * n;
}

void initFact()
{
	for(int i= 0; i <= MAX_N; i++)
		fact[i] = factorial(i);
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
		sum += powInt(x, i) / fact[i];
		
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
			= absErrMatrix.val[i][j]*fact[i +1]
			/powInt( absErrMatrix.val[0][j], i + 1 );
	
	printDatMatD( absErrMatrix, "e_approx.dat",   "%.14lf ", true);
	printDatMatD( relErrValues, "rel_errors.dat", "%.14lf ", true);
	
	freeAllMatD(absErrMatrix, relErrValues, NULL_MAT);
}


int main()
{
	initFact();
	analysis(MAX_N, N_INTERVALS);
	
	return 0;
}
