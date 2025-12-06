// plot.c

#define _DEFAULT_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include "comp_physics.h"
#include <stdbool.h>
#include <string.h>

FILE * gp_open(  )
{
	FILE * gp;
	if ( !(gp = popen("gnuplot", "w")) )
		raiseErr("Cannot run gnuplot");
	return gp;
}

void gp_term_def( FILE * gp, char * title )
{
	fprintf( gp, "set terminal qt\n" );
	if (title)
		fprintf( gp, "set title '%s'\n", title );
}

void gp_end( FILE * gp )
{
	fprintf( gp, "pause mouse close\n" );
	fflush(gp);
	pclose(gp);
}

void gp_prt_carr( FILE * gp, double * xcarr, double * ycarr, int length )
{
	for( int j = 0; j < length; j++ )
		fprintf( gp, "%g %g\n", xcarr[j], ycarr[j] );
	fprintf( gp, "e\n" );
}

static inline bool is_requested_present( char ** s_arr, int num )
{
	if( s_arr )
		if( s_arr[num] )
			return true;
	
	return false;
}

void gp_set_plot( FILE * gp, int num, char ** labels, char ** styles )
{
	fprintf(gp, "plot ");
    for(int i = 0; i < num; i++)
    {
		
		fprintf(gp, "'-' with %s title '%s'%s", 
			is_requested_present( styles, i ) ? styles[i]: "lines",
			is_requested_present( labels, i ) ? labels[i]: "" ,
			( i < num -1 ) ? ", ": "\n" );
	}
}

void gp_set_logscale( FILE * gp, bool x, bool y )
{
	if( x )
		fprintf( gp, "set logscale x\n" );
	if( y )
		fprintf( gp, "set logscale y\n" );
}

void gp_axes_labels( FILE * gp, char * x_label, char * y_label )
{
	if( x_label )
		fprintf( gp, "set xlabel '%s'\n", x_label );
		
	if( y_label )
		fprintf( gp, "set ylabel '%s'\n", y_label );
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
	
	FILE * gp = gp_open();
		
	gp_term_def( gp, NULL );
	fprintf(gp, "set grid nopolar\n");
	gp_set_plot( gp, 1, NULL, NULL );
	gp_prt_carr( gp, x->val, y->val, x->length );
	
	gp_end(gp);
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
	
	gp_term_def( gp, plot_info.title );
		
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
		gp_prt_carr( gp, x_axis->val, y_values[i].val, x_axis->length );
		
	gp_end( gp );
}

static inline void append_and_advance( char **dest, char c ){ **dest = c; ++*dest; }

static inline bool is_c_fine( char c, char stop )
{
	return c != stop && c != EOF;
}

static inline bool read_until_char( FILE * file, char stop, char ** dest )
{
	char c;
	while( is_c_fine( c = (char)getc(file), stop ) )
		append_and_advance( dest, c );
		
	if (c == EOF)
		return true;
		
	return false;
}

static inline bool parse_line( FILE * cfg, char ** dest )
{

	if( read_until_char( cfg, '=', dest ) )
		return false;
	
	append_and_advance( dest, '=' );
	append_and_advance( dest, '\'' );
			
	read_until_char( cfg, '\n', dest );
	
	append_and_advance(dest, '\'');
	append_and_advance(dest, ';');
	append_and_advance(dest, ' ');
	
	return true;
} 

static inline void config_parser( const char * filename, char * dest )
{
	FILE * config = openFile( filename, "r" );
	
	while( parse_line( config, &dest) );
	
	*dest = '\0';
	
	fclose(config);
}

static inline bool are_same_length(const VectorD * vecs[], int num_vecs)
{
	for( int i = 1; i < num_vecs; i++ )
		if( vecs[i]->length != vecs[0]->length )
			return false;
	return true;
}

static inline void vecs_to_file( FILE * file, const VectorD * vecs[], int num_vecs )
{
	if( !are_same_length(vecs, num_vecs) )
		raiseErr("vecs are not all same length.");
	
	const int vec_length = vecs[0]->length;
	
	for( int i = 0; i < vec_length; i++ )
	{
		for( int j = 0; j < num_vecs; j++ )
		{
			fprintf( file, "%16g ", vecs[j]->val[i]);
		}
		putc('\n', file);
	}
}

void two_carr_to_file( const char * filename, const VectorD * vec1, const VectorD * vec2 )
{
	FILE * data = openFile( filename, "w" );
	const VectorD * vecs[2] = { vec1, vec2 };
	vecs_to_file( data, vecs, 2 );
	fclose(data);
}

void tmplot_2vecs( const char * cfg_name, const VectorD * xaxis, const VectorD * yaxis)
{
	char command [2000];
	char parsed_cfg[1000];
	
	char cfg_path[100];
	sprintf(cfg_path, "conf/%s.cfg", cfg_name);
	
	char data_path[100];
	sprintf( data_path, "data/%s.dat", cfg_name);

	const char * tmp_path = "../../gp_templates";

	config_parser( cfg_path, parsed_cfg );
	
	two_carr_to_file( data_path, xaxis, yaxis );
	
	sprintf( command,"gnuplot -e \"%sMY_DATAFILE=\'%s\'; MODE=\'one\'; "
			"TMP_PATH=\'%s\'\" ../../gp_templates/main.gp", 
			parsed_cfg, data_path, tmp_path );
//	printf(command);
	system(command);
}



