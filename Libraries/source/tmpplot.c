// tmpplot.c

#include <stdlib.h>
#include <stdio.h>
#include "comp_physics.h"
#include <stdbool.h>
#include <string.h>

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

static inline void carrs_to_file( FILE * file, const VectorD * vecs[], int num_vecs )
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

void tmplot_2carr( const char *cfg_name, double *xaxis, double *yaxis, int length)
{
	char command [2000];
	char parsed_cfg[1000];
	
	char cfg_path[100];
	sprintf(cfg_path, "conf/%s.cfg", cfg_name);
	
	char data_path[100];
	sprintf( data_path, "data/%s.dat", cfg_name);

	const char *tmp_path = "../../gp_templates";

	config_parser( cfg_path, parsed_cfg );
	
	two_carr_to_file( data_path, xaxis, yaxis );
	
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
	
	two_carr_to_file( data_path, xaxis, yaxis );
	
	sprintf( command,"gnuplot -e \"%sMY_DATAFILE=\'%s\'; MODE=\'one\'; "
			"TMP_PATH=\'%s\'\" ../../gp_templates/main.gp", 
			parsed_cfg, data_path, tmp_path );
//	printf(command);
	system(command);
}













