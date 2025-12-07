#include "comp_physics.h"
#include <math.h>

#define MAX_DEGREE_TEST 10

double leg2roots[2];
double leg3roots[3];
double leg4roots[4];
double leg5roots[5];
double *leg_roots[]= {leg2roots, leg3roots, leg4roots, leg5roots};

void setup_leg_roots(){
	leg2roots[0]=-sqrt(3)/3;
	leg2roots[1]= sqrt(3)/3;
	leg3roots[0]= sqrt(15)/5;
	leg3roots[2]=-sqrt(15)/5;
	leg3roots[1]=0;
	leg4roots[0]=-sqrt( 525 + 70 * sqrt(30) )/35;
	leg4roots[3]= sqrt( 525 + 70 * sqrt(30) )/35;
	leg4roots[1]=-sqrt( 525 - 70 * sqrt(30) )/35;
	leg4roots[2]= sqrt( 525 - 70 * sqrt(30) )/35;

}

int main()
{
	setup_leg_roots();

	
	for( int deg = 1; deg <= MAX_DEGREE_TEST; deg++ )
	{
		printf("n = %d\n", deg);
		for( int root = 1; root <= deg; root++ )
		{
			printf( "k = %2d: %16g\n", root, legendre_root(deg, root) );
		}
		printf("-------------------------------------\n");
	}

	
	exit(EXIT_SUCCESS);
}

