#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "explot.h"
#include "genutil.h"

static int gb_counter=0;

void eplot_set_common(FILE *gp)
{
	const char * s = "set encoding utf8\n"
		"set mxtics 1\n"
		"set mytics 1\n"
		"set grid xtics ytics mxtics mytics lw 0.4 lc rgb '#AAAAAA'\n"
		"set border 3 lw 0.5 lc rgb \"black\"\n"
		"set xtics nomirror scale 1\n"
		"set ytics nomirror scale 1\n"
		"set key box lw 0.5 lc rgb \'black\' textcolor rgb \'black\' "
		"spacing 1.5 opaque width 0.5\n"
		// global line presets	
		"set linetype 1 lc rgb \"#66C2A5\"\n"
		"set linetype 2 lc rgb \"#FC8D62\"\n"
		"set linetype 3 lc rgb \"#8DA0CB\"\n"
		"set linetype 4 lc rgb \"#E78AC3\"\n"
		"set linetype 5 lc rgb \"#A6D854\"\n"
		"set linetype 6 lc rgb \"#FFD92F\"\n"
		"set linetype 7 lc rgb \"#E5C494\"\n"
		"set linetype 8 lc rgb \"#B3B3B3\"\n";
	if(fputs(s, gp)==EOF)
		raiseErr("could not write common settings");
}

//lc rgb '#66C2A5'

Global null_global_settings()
{
	Global settings = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
			nullptr, nullptr};
	return settings;
}

DatasetDesc null_datadesc()
{
	DatasetDesc desc={nullptr, nullptr, nullptr};
	return desc;
}

void eplot_set_global(FILE *gp, const Global *settings )
{
	if(settings->title)
		fprintf(gp, "set title \"%s\" font \',24\'\n", settings->title);
	
	if(settings->xlabel)
		fprintf(gp, "set xlabel \"%s\" font \',20\'\n", settings->xlabel);
	else
		fprintf(gp, "set xlabel \"x\" font \',20\'\n");
		
	if(settings->ylabel)
		fprintf(gp, "set ylabel \"%s\" font \',20\'\n", settings->ylabel);
	else
		fprintf(gp, "set ylabel \"y\" font \',20\'\n");
	
	if(settings->logscale)
		fprintf(gp, "set logscale %s\n", settings->logscale);
	
	if(settings->key)
		fprintf(gp, "set key %s\n", settings->key);
		
	if(settings->additional_set)
		fprintf(gp, "%s\n", settings->additional_set);
}

void eplot_set_terminal_and_plot(FILE *gp, char *plot_command, const char *output_name)
{
	char qt[5'000];
	sprintf(qt,
		"set terminal qt size 1020,700 font \'CMU Serif, 16\' enhanced\n"
		"set for [i=1:8] linetype i lw 4 pt 7 ps 1.5\n"
		"%s\n"
		"pause mouse close\n", plot_command);
	fputs(qt, gp);
	
	if(output_name!=nullptr)
	{
		char pdfcairo[5'000];
		char output_path[100];
		sprintf(output_path, "\'plot/%s.pdf\'", output_name);
		sprintf(pdfcairo,
			"set terminal pdfcairo size 6.4,4 font \"CMU Serif,16\" enhanced\n"
			"set for [i=1:8] linetype i lw 4 pt 7 ps 1.5\n"
			"set output %s\n"
			"%s\n"
			"unset output\n", output_path, plot_command );
		fputs(pdfcairo, gp);
	}
}

static inline void carrs_to_file( FILE *file, double **carrs, int num_carrs,
									int length )
{	
	for( int i = 0; i < length; i++ )
	{
		for( int j = 0; j < num_carrs; j++ )
		{
			fprintf( file, "%.16g ", carrs[j][i]);
		}
		putc('\n', file);
	}
}

static inline void carrs_to_filename( const char *filename, double **carrs,
										int length, int num_carrs)
{
	FILE *data = openFile( filename, "w" );
	carrs_to_file( data, carrs, num_carrs, length );
	fclose(data);
}

static inline void plot_command_helper_2carr(char *plot_command, 
		const DatasetDesc *data_settings, char *data_filename)
{
	char *label=(data_settings->label) ? data_settings->label : "";
	char *style=(data_settings->style) ? data_settings->style : "lp";
	
	
	sprintf(plot_command, "plot '%s' with %s ls 1 title '%s'\n", data_filename,
				style, label);
}

static inline char *get_build_name(const Global *gb_settings)
{return (gb_settings->build_name) ? gb_settings->build_name :"my_plot";}

static inline FILE *set_and_open_script(char *build_name, char *script_path){
	sprintf(script_path, "script/%s.gp", build_name);
	FILE *gp=openFile(script_path, "w");
	return gp;
}

void eplot_2carr(double *xaxis, double *yaxis, int size, const Global *gb_settings,
					const DatasetDesc *data_settings)
{
	char *build_name=get_build_name(gb_settings);
	
	char script_path[100]; 
	FILE *gp=set_and_open_script(build_name, script_path);
	
	eplot_set_common(gp);
	
	eplot_set_global(gp, gb_settings);

	char data_path[100];
	sprintf(data_path, "data/%s_%d_0.dat", build_name, gb_counter);
	double *carrs[2]={xaxis,yaxis};
	carrs_to_filename(data_path, carrs, size, 2);

	char plot_command[2'000];
	plot_command_helper_2carr(plot_command, data_settings, data_path);
//	fputs(plot_command, gp);

	eprint("%s", plot_command);
	eplot_set_terminal_and_plot(gp, plot_command, gb_settings->output_name);

	fflush(gp);
	fclose(gp);

	char shell_command[200];
	sprintf(shell_command, "gnuplot %s", script_path);
	system(shell_command);
	gb_counter++;
}



void multi_data_plt_command_helper(const SeriesSpec *series, char *plot_command,
										const char *build_name, int i)
{
	char data_path[100];
	sprintf(data_path, "data/%s_%d_%d.dat", build_name, gb_counter, i);
	double *carrs[2]={series->x,series->y};
	carrs_to_filename(data_path, carrs, series->size, 2);
	
	char *label=(series->label) ? series->label : "";
	char *style=(series->style) ? series->style : "lp";
	
	char plot_line[400];
	sprintf(plot_line, "\"%s\" with %s title \"%s\"", data_path,
				style, label);
	strcat(plot_command, plot_line);
}

void eplot_multi(const SeriesSpec *series, int count, const Global *gb_settings)
{
	char *build_name=get_build_name(gb_settings);
	
	char script_path[100]; 
	FILE *gp=set_and_open_script(build_name, script_path);
	
	eplot_set_common(gp);
	
	eplot_set_global(gp, gb_settings);
	
	char plot_command[2'000]={0};
	strcat(plot_command, "plot ");
	
	for(int i=0; i<count; i++)
	{	
		multi_data_plt_command_helper(series+i, plot_command, build_name, i);
		strcat(plot_command, (i==count-1)?"\n":", \\" "\n");
	}
	
	eprint("%s", plot_command);
	eplot_set_terminal_and_plot(gp, plot_command, gb_settings->output_name);

	fflush(gp);
	fclose(gp);

	char shell_command[200];
	sprintf(shell_command, "gnuplot %s", script_path);
	system(shell_command);
	
	gb_counter++;
}


