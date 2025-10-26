#include <stdlib.h>
#include <stdio.h>
#include <DoubleUtilities.c>
#include <LinearSystems.c>
#include <math.h>

#define PI 3.14159265358979323846

double logarithm(double x){return log(x);}

double cst_1( double x )
{
	UNUSED(x);
	return 1;
}

void kepler_law_int()
{
	ArrayDouble distances = buildArrD( 8, 
		57.59,
		108.11,
		149.57,
		227.84,
		778.14,
		1427.0,
		2870.3,
		4499.9);
		
	ArrayDouble periods = buildArrD( 8,
		87.99,
		224.7,
		365.26,
		686.98,
		4332.4,
		10759.0,
		30684.0,
		60188.0);
		
	ArrayDouble log_periods = allocArrD( 8 );
	for( int i = 0; i < 8; i++)
		log_periods.val[i] = log(periods.val[i]);
	
	Linear_Model model = alloc_linear_model( 2 );
	model.funcs[0] = cst_1;
	model.funcs[1] = logarithm;
	
	ArrayDouble result =  linear_least_square_fitting( distances, log_periods,
		model);
	
	printArrDPar( result, "%lf, " );
	
	free(model.funcs);
	freeAllArrD( distances, periods, result, NULL_ARR );
	
	// FARE PLOT
}

int main()
{
	kepler_law_int();
	
	return 0;
}
