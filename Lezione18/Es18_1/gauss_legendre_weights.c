#include "comp_physics.h"

#define MAX_DEGREE_TEST 10

int main()
{
	
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

