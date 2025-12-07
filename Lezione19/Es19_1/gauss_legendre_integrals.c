#include "comp_physics.h"
#include <math.h>

double f1( double x, void * unused ){
	UNUSED(unused);
	return exp(-4*x);}
	
Par_Func f1_p = { f1, 0, 0 };

int main()
{
	for( int n = 4; n<=40; n+=4 )
		fprintf(stderr, "%16g\n%16g\n", int_gauss_legendre( &f1_p, n ), sinh( 4 )/2 );
	
	return 0;
}
