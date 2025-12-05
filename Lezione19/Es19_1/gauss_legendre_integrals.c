#include "comp_physics.h"
#include <math.h>

double f1( double x, void * unused ){
	UNUSED(unused);
	return exp(-4*x);}
	
Par_Func f1_p = { f1, 0, 0 };

int main()
{
	printf( "%g\n%g\n", int_gauss_legendre( &f1_p, 40 ), sinh( 4 )/2 );
	
	return 0;
}
