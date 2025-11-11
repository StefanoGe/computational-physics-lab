// plot.c

#define _DEFAULT_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include "comp_physics.h"

static inline FILE * gp_open(  )
{
	FILE * gnuplot;
	if ( !(gnuplot = popen("gnuplot", "w")) )
		raiseErr("Cannot run gnuplot");
	return gnuplot;
}

PlotInfo plot_info_init(int num)
{
	PlotInfo plot_info;
	plot_info.title = NULL;
	plot_info.labels = malloc( num * sizeof(char *) );
	for(int i = 0; i < num; i++)
		plot_info.labels[i] = NULL;
		
	plot_info.num = num;		
	return plot_info;
}

void plot_2vecs( const VectorD * x, const VectorD * y )
{
	if( x->length != y->length )
		raiseErr( "x and y must be of same length" );
	
	FILE * gnuplot = gp_open();
		
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

void plot_mult_vecs_var( int num, const VectorD * x_axis, ... )
{
	va_list list;
	va_start(list, x_axis);
	
	FILE * gp = gp_open();
	
	fprintf( gp, "set terminal qt\n" );
	fprintf(gp, "plot ");
	
    for(int i = 0; i < num; i++)
    {
        fprintf(gp, "'-'%s", ( i < num -1 ) ? ", ": "\n" );
	}
	
	VectorD * curr_vec;
	
	for( int i = 0; i < num; i++ )
	{
		if( !( curr_vec = va_arg( list, VectorD * ) ) )
			raiseErr( "cannot get Vector * from list" );
		for( int j = 0; j < curr_vec->length; j++ )
			fprintf( gp, "%g %g\n", x_axis->val[j], curr_vec->val[j] );
		fprintf( gp, "e\n" );
	}
	va_end(list);
	
	fflush( gp );
	pclose( gp );
}

void plot_mult_vecs( const VectorD * x_axis, const VectorD y_values[], 
					int num_vec, PlotInfo plot_info )
{
	FILE * gp = gp_open();
	
	fprintf( gp, "set terminal qt\n" );

	if( plot_info.title )
		fprintf( gp, "set title %s\n", plot_info.title );
		
	fprintf( gp, "set key box\n" );
	
	fprintf(gp, "plot ");
	
    for(int i = 0; i < num_vec; i++)
    {
		if( plot_info.labels )
			fprintf(gp, "'-' title '%s'%s", 
				plot_info.labels[i] ? plot_info.labels[i]: "" ,
				( i < num_vec -1 ) ? ", ": "\n" );
		else
			fprintf(gp, "'-'%s", ( i < num_vec -1 ) ? ", ": "\n" );
	}
	
	for( int i = 0; i < num_vec; i++ )
	{
		for( int j = 0; j < y_values[i].length; j++ )
			fprintf( gp, "%g %g\n", x_axis->val[j], y_values[i].val[j] );
		fprintf( gp, "e\n" );
	}
	fprintf( gp, "pause mouse close\n" );
	
	fflush( gp );
	pclose( gp );
}

