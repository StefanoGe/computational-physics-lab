// Prova

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comp_physics.h"

#define MAX_N 4
#define N_INTERVALS 10
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

double naiveExp( double x, int N )
{
	double sum = 0;
	
	for(int i=0; i<=N; i++)
		sum += pown(x, i) / factorial(i);
		
	return sum;
}

double smartExp( double x, int N )
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
	return fabs( smartExp( x, N ) - exp( x ) );
}

void analysis ( int max_N, int n_intervals )
{
	const int n_xvalues = n_intervals + 1;
	ArrayDouble xAxis = linspaceD( 0, 1, n_xvalues);
	MatrixDouble errs = allocMatD(max_N, n_xvalues);
	
	// We compute the relative error matrix.
	for( int i = 1; i <= max_N; i++ )
		for(int j = 1; j < n_xvalues; j++)
		{
			errs.val[i-1][j] = absoluteError( xAxis.val[j], i);
			printf("x = %lf - exp = %lf - err = %lf\n", 
				xAxis.val[j], smartExp(xAxis.val[j], i), errs.val[i-1][j]);
			
		}

	// Plot
	tmplot_carrs( "plot1", xAxis.val, errs.val, n_xvalues, max_N );
	freeArrD(xAxis);
	freeAllMatD(errs, NULL_MAT);
}


int main()
{
	analysis(MAX_N, N_INTERVALS);
	
	return 0;
}
