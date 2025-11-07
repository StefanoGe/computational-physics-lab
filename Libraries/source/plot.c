// plot.c

#include "comp_physics.h"
#include <stdlib.h>

void plot_2vecs( const VectorD * x, const VectorD * y )
{
	if( x->length != y->length )
		raiseErr( "x and y must be of same length" );
	
	FILE *gnuplot = popen("gnuplot", "w");
	
	fprintf(gnuplot, "plot '-'\n");
	for (int i = 0; i < x->length; i++)
		fprintf(gnuplot, "%g %g\n", x->val[i], y->val[i]);
	fprintf(gnuplot, "e\n");
	fflush(gnuplot);
}

void plot_func( VectorD * domain, Func_Ptr func )
{
	const int length = domain->length;
	VectorD y = init_vec_length( length );
	for( int i = 0; i < length; i ++ )
		y.val[ i ] = func( domain ->val[i] );
		
	eprint("fine");
	
	plot_2vecs( domain, &y );
	
	free_vecD( &y );
}
