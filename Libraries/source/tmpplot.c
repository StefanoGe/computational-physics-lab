// tmpplot.c

#include <stdlib.h>
#include <stdio.h>
#include "comp_physics.h"
#include <stdbool.h>
#include <string.h>
#include <ctype.h>

#define MAX_SIZE_COMMAND 11'000
#define MAX_SIZE_CONFIG 10'000
const char *tmp_path = "../../gp_templates";

static inline bool read_line(FILE *file, char *buf, int cap) 
{
    if (!fgets(buf, cap, file)) return false;
    return true;
}

static bool is_comment_or_empty(const char *line) 
{
	// Skips starting spaces or newline
    while (*line && isspace((int)*line)) line++;
    return (*line == '\0' || *line == '#');
}

static bool parse_key_value(const char *line, char *dest, int dest_cap) {
    const char *eq = strchr(line, '=');
    if (!eq) return false;

    // Extract key
    char key[256];
    size_t klen = eq - line;
    if (klen >= sizeof key) return false;

    memcpy(key, line, klen);
    key[klen] = '\0';

    // Trim trailing spaces in key
    while (klen > 0 && isspace((unsigned char)key[klen - 1])) {
        key[--klen] = '\0';
    }

    // Extract value (raw)
    const char *value = eq + 1;
    while (*value && isspace((unsigned char)*value)) value++;

    // Trim trailing newline/spaces
    size_t vlen = strlen(value);
    while (vlen > 0 && isspace((unsigned char)value[vlen - 1])) vlen--;

    // Append "KEY='VALUE'; "
    int wrote = snprintf(dest, dest_cap, "%s='%.*s'; ", key, (int)vlen, value);
    return (wrote >= 0 && (int)wrote < dest_cap);
}

void config_parser(const char *filename, char *out, size_t out_cap) {
    FILE *file = openFile(filename, "r");

    char line[1024];
    size_t used = 0;

    while (read_line(file, line, sizeof line)) {
        if (is_comment_or_empty(line)) continue;

        // Build into a temp buffer, then append
        char tmp[1024];
        if (parse_key_value(line, tmp, sizeof tmp)) 
        {
            size_t len = strlen(tmp);
            if (used + len + 1 < out_cap) 
            {
                memcpy(out + used, tmp, len);
                used += len;
            }
        }
    }

    out[used] = '\0';
    fclose(file);
}
/*
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

static inline void config_parser_old( const char * filename, char * dest )
{
	FILE * config = openFile( filename, "r" );
	
	while( parse_line( config, &dest) );
	
	*dest = '\0';
	
	fclose(config);
}
*/
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
			fprintf( file, "%16g ", carrs[j][i]);
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

static void command_build_helper(const char *cfg_name, 
				char *cfg_path, char *data_path, char *parsed_cfg)
{
	sprintf(cfg_path, "conf/%s.cfg", cfg_name);
	sprintf( data_path, "data/%s.dat", cfg_name);

	config_parser( cfg_path, parsed_cfg, MAX_SIZE_CONFIG);
}

static void plot_command( const char *parsed_cfg,
	const char *data_path, const char *tmp_path, const char *mode)
{
	char command [MAX_SIZE_COMMAND];
	sprintf( command,"gnuplot -e \"%sMY_DATAFILE=\'%s\'; MODE=\'%s\'; "
			"TMP_PATH=\'%s\'\" ../../gp_templates/main.gp", 
			parsed_cfg, data_path, mode, tmp_path );
//	printf(command);
	system(command);
}

void tmplot_2carr( const char *cfg_name, const double *xaxis, 
							const double *yaxis, int length)
{
	char parsed_cfg[MAX_SIZE_CONFIG];
	char cfg_path[100];
	char data_path[100];
	
	command_build_helper( cfg_name, cfg_path, data_path, parsed_cfg );
	
	two_carr_to_filename( data_path, xaxis, yaxis, length );
	
	plot_command( parsed_cfg, data_path, tmp_path, "one" );
}

void tmplot_2vecs( const char * cfg_name, const VectorD * xaxis, const VectorD * yaxis)
{
	char parsed_cfg[MAX_SIZE_CONFIG];
	char cfg_path[100];
	char data_path[100];
	
	command_build_helper( cfg_name, cfg_path, data_path, parsed_cfg );

	two_vecs_to_filename( data_path, xaxis, yaxis );

	plot_command( parsed_cfg, data_path, tmp_path, "one" );
}

void tmplot_carrs( const char *cfg_name, const double *xaxis, 
							double **yaxes, int length, int num_yaxes)
{
	char parsed_cfg[MAX_SIZE_CONFIG];
	char cfg_path[100];
	char data_path[100];
	
	command_build_helper( cfg_name, cfg_path, data_path, parsed_cfg );
	
	const int num_carrs = num_yaxes+1;
	const double **all_carrs = malloc( sizeof(double*) * num_carrs);
	all_carrs[0]=xaxis;
	for(int i = 1; i<num_carrs; i++)
		all_carrs[i] = yaxes[i-1];
	
	carrs_to_filename( data_path, all_carrs, length, num_carrs );
	free(all_carrs);
	
	plot_command( parsed_cfg, data_path, tmp_path, "multi" );
}











