#include <stdlib.h>
#include <stdio.h>
#include "comp_physics.h"
#include <math.h>

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

int main()
{
	test_min_x();
	
	return 0;
}





