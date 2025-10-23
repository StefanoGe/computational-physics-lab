#include <stdlib.h>
#include <stdio.h>
#include <DoubleUtilities.c>
#include <LinearSystems.c>
#include <math.h>

#define PI 3.14159265358979323846


#define MONO_LIST \
    X(0) \
    X(1) \
    X(2) \
    X(3) \
    X(4)

#define X(N) double mono_##N(double x) { return pow(x, N); }
MONO_LIST
#undef X

double logarithm( double x){return log(x);}

void test_min_x(  )
{
	MatrixDouble test_mat = readMatD( "data/test_mat_ls.txt", 3, 2 );
	ArrayDouble test_b = buildArrD( 3, 1.0, -5.0, 6.0 );
	ArrayDouble minimum_x = min_x( test_mat, test_b );
	
	printf("x which minimizes required least square system:\n");
	printArrDPar( minimum_x, "%lf " );
	
	freeAllArrD( test_b, minimum_x, NULL_ARR );
	freeMatD( test_mat );
	
	
	// right result when compared to result from website.
}

double ret_1( double x )
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
	model.funcs[0] = ret_1;
	model.funcs[1] = logarithm;
	
	ArrayDouble result =  linear_least_square_fitting( distances, log_periods,
		model);
	
	printArrDPar( result, "%lf, " );
	
	free(model.funcs);
	freeAllArrD( distances, periods, result, NULL_ARR );
	
	// FARE PLOT
}

double e_sin(double t)
{
	return exp(sin(t-1));
}

void approx_e_sin()
{
	//ArrayDouble x_data = linspaceD( 2 * PI / 60, 1, 60 );
	//ArrayDouble y_data = build_from( x_data, e_sin );
	
	
	printf("%lf\n", mono_3(2));
	
}

int main()
{
	//eprint("%lf\n", log(60188));
	//test_func_point();
	// test_min_x();
	// kepler_law_int();
	
	approx_e_sin();
	
	return 0;
}





