// plot.c

#define _DEFAULT_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include "comp_physics.h"


void plot_2vecs( const VectorD * x, const VectorD * y )
{
	if( x->length != y->length )
		raiseErr( "x and y must be of same length" );
	
	FILE *gnuplot;
	if ( !(gnuplot = popen("gnuplot", "w")) )
		raiseErr("Cannot run gnuplot");
		
	fprintf( gnuplot, "set terminal qt\n" );
	
	fprintf(gnuplot, "plot '-'\n");
	for (int i = 0; i < x->length; i++)
		fprintf(gnuplot, "%g %g\n", x->val[i], y->val[i]);
	fprintf(gnuplot, "e\n");
	
	fprintf( gnuplot, "pause mouse close\n" );

	fflush(gnuplot);
	pclose(gnuplot);
}

void plot_func( VectorD * domain, Func_Ptr func )
{
	const int length = domain->length;
	VectorD y = init_vec_length( length );
	for( int i = 0; i < length; i ++ )
		y.val[ i ] = func( domain ->val[i] );
			
	plot_2vecs( domain, &y );
	
	free_vecD( &y );
}
