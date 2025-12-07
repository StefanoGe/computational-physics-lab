// tmpplot.c

#include <stdlib.h>
#include <stdio.h>
#include "comp_physics.h"
#include <stdbool.h>
#include <string.h>

#define ALL_GOOD 1
#define REACHED_EOF 2
#define FOUND_NEWLINE 3
#define FOUND_HASH 4

static inline void append_and_advance( char **dest, char c ){ **dest = c; ++*dest; }

static inline bool is_c_fine( char c, char stop, int *state )
{
	return c != stop && c != EOF;
}

static inline bool read_until_char( FILE * file, char stop, char ** dest )
{
	char c;
	int state=ALL_GOOD;
	if( c == '\n' )
		return true;
	if( c == '#')
	{
		while( c!= '\n' || c != EOF)
		{
			
		}
	}
	
	while( is_c_fine( c = (char)getc(file), stop, &state ) )
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

static inline void carrs_to_file( FILE *file, const double **carrs, int num_carrs,
									int length )
{	
	for( int i = 0; i < length; i++ )
	{
		for( int j = 0; j < num_carrs; j++ )
		{
			fprintf( file, "%16g ", carrs[i][j]);
		}
		putc('\n', file);
	}
}

static inline void carrs_to_filename( const char *filename, const double **carrs,
										int length, int num_carrs)
{
	FILE * data = openFile( filename, "w" );
	carrs_to_file( data, carrs, num_carrs, length );
	fclose(data);
}

static inline void two_vecs_to_filename( const char *filename, const VectorD *vec1, const VectorD *vec2 )
{
	if( vec1->length != vec2->length )
		raiseErr("vecs are not all same length.");
	FILE * data = openFile( filename, "w" );
	const double * vecs[2] = { vec1->val, vec2->val };
	carrs_to_file( data, vecs, 2, vec1->length );
	fclose(data);
}

static inline void two_carr_to_filename( const char *filename, const double *carr1, 
					const double *carr2, int length )
{
	FILE * data = openFile( filename, "w" );
	const double *carrs[2] = { carr1, carr2 };
	carrs_to_file( data, carrs, 2, length );
	fclose(data);
}

void tmplot_2carr( const char *cfg_name, const double *xaxis, 
							const double *yaxis, int length)
{
	char command [2000];
	char parsed_cfg[1000];
	
	char cfg_path[100];
	sprintf(cfg_path, "conf/%s.cfg", cfg_name);
	
	char data_path[100];
	sprintf( data_path, "data/%s.dat", cfg_name);

	const char *tmp_path = "../../gp_templates";

	config_parser( cfg_path, parsed_cfg );
	
	two_carr_to_filename( data_path, xaxis, yaxis, length );
	
	sprintf( command,"gnuplot -e \"%sMY_DATAFILE=\'%s\'; MODE=\'one\'; "
			"TMP_PATH=\'%s\'\" ../../gp_templates/main.gp", 
			parsed_cfg, data_path, tmp_path );
//	printf(command);
	system(command);
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
	
	two_vecs_to_filename( data_path, xaxis, yaxis );
	
	sprintf( command,"gnuplot -e \"%sMY_DATAFILE=\'%s\'; MODE=\'one\'; "
			"TMP_PATH=\'%s\'\" ../../gp_templates/main.gp", 
			parsed_cfg, data_path, tmp_path );
//	printf(command);
	system(command);
}

void tmplot_carrs( const char *cfg_name, const double *xaxis, 
							const double **yaxes, int length, int num_yaxes)
{
	char command [2000];
	char parsed_cfg[1000];
	
	char cfg_path[100];
	sprintf(cfg_path, "conf/%s.cfg", cfg_name);
	
	char data_path[100];
	sprintf( data_path, "data/%s.dat", cfg_name);

	const char *tmp_path = "../../gp_templates";

	config_parser( cfg_path, parsed_cfg );
	
	const int num_carrs = num_yaxes+1;
	const double **all_carrs = malloc( sizeof(double*) * num_carrs);
	all_carrs[0]=xaxis;
	for(int i = 1; i<=num_carrs; i++)
		all_carrs[i] = yaxes[i-1];
	
	carrs_to_filename( data_path, all_carrs, length, num_carrs );
	free(all_carrs);
	
	sprintf( command,"gnuplot -e \"%sMY_DATAFILE=\'%s\'; MODE=\'multi\'; "
			"TMP_PATH=\'%s\'\" ../../gp_templates/main.gp", 
			parsed_cfg, data_path, tmp_path );
//	printf(command);
	system(command);
}











