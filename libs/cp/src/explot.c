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
	snprintf(qt, sizeof qt,
		"set terminal qt size 1020,700 font \'CMU Serif, 16\' enhanced\n"
		"set for [i=1:8] linetype i lw 4 pt 7 ps 1.6\n"
		"%s\n"
		"pause mouse close\n", plot_command);
	fputs(qt, gp);
	
	if(output_name!=nullptr)
	{
		char pdfcairo[5'000];
		char output_path[100];
		sprintf(output_path, "\'plot/%s_%d.pdf\'", output_name, gb_counter);
		sprintf(pdfcairo,
			"set terminal pdfcairo size 6.4,4 font \"CMU Serif,16\" enhanced\n"
			"set for [i=1:8] linetype i lw 4.5 pt 7 ps 0.65\n"
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
	sprintf(script_path, "script/%s_.gp", build_name);
	FILE *gp=openFile(script_path, "w");
	return gp;
}

void eplot_2carr(double *xaxis, double *yaxis, int size, const Global *gb_settings,
					const DatasetDesc *data_settings)
{
	char build_name[100];
	if(gb_settings->build_name)
		snprintf(build_name, sizeof build_name, "%s_%d", gb_settings->build_name,
		gb_counter);
	else
		snprintf(build_name, sizeof build_name, "my_plot_%d", gb_counter);

	
	char script_path[100]; 
	FILE *gp=set_and_open_script(build_name, script_path);
	
	eplot_set_common(gp);
	
	eplot_set_global(gp, gb_settings);

	char data_path[200];
	sprintf(data_path, "data/%s_0.dat", build_name);
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
										const char *build_name, int i,size_t *len,
										size_t plot_command_max_size)
{
	char data_path[100];
	snprintf(data_path,sizeof data_path, "data/%s_%d.dat", build_name, i);
	double *carrs[2]={series->x,series->y};
	carrs_to_filename(data_path, carrs, series->size, 2);
	
	char *label=(series->label) ? series->label : "";
	char *style=(series->style) ? series->style : "lp";
	char color[50]="";
	if(series->color)
		snprintf(color, sizeof color, " lc %s", series->color);
	
	char plot_line[400];
	snprintf(plot_line, sizeof plot_line,"\"%s\" with %s title \"%s\"%s", 
				data_path,
				style, label, color);
	*len+=snprintf(plot_command+*len,plot_command_max_size-*len, plot_line);
}

void eplot_multi(const SeriesSpec *series, int count, 
	const Global *gb_settings)
{
	eplot_multi_options(series, count, gb_settings, true);
}

void eplot_multi_options(const SeriesSpec *series, int count, 
	const Global *gb_settings, bool print_plot_command)
{
	char build_name[200];
	if(gb_settings->build_name)
		snprintf(build_name, sizeof build_name, "%s_%d", gb_settings->build_name,
		gb_counter);
	else
		snprintf(build_name, sizeof build_name, "my_plot_%d", gb_counter);
	
	char script_path[100]; 
	FILE *gp=set_and_open_script(build_name, script_path);
	
	eplot_set_common(gp);
	
	eplot_set_global(gp, gb_settings);
	
	char plot_command[2'000];
	size_t len=snprintf(plot_command,sizeof plot_command, "plot ");
	
	for(int i=0; i<count; i++)
	{	
		multi_data_plt_command_helper(series+i, plot_command, build_name, i,&len,
			sizeof plot_command);
		len+=snprintf(plot_command+len,sizeof plot_command - len,
												(i==count-1)?"\n":", \\" "\n");
	}
	if(print_plot_command)
		eprint("%s", plot_command);
	eplot_set_terminal_and_plot(gp, plot_command, gb_settings->output_name);

	fflush(gp);
	fclose(gp);

	char shell_command[200];
	snprintf(shell_command, sizeof shell_command, "gnuplot %s", script_path);
	system(shell_command);
	
	gb_counter++;
}

void anim_data_plt_command_helper1(char *plot_command,
										const char *build_name, int i,
										size_t *len,
										size_t plot_command_max_size,
										bool is_static)
{
	*len+=snprintf(plot_command+*len,plot_command_max_size-*len,
				"file%d = sprintf(\"data/%s_%d_%d_",
				i, build_name,gb_counter,i);
	*len+=snprintf(plot_command+*len,plot_command_max_size-*len,"%%");
	*len+=snprintf(plot_command+*len,plot_command_max_size-*len,
				(is_static?"05d.dat\", 0)\n":"05d.dat\", i)\n"));
}

void anim_data_plt_command_helper2(const AnimSpec *series, char *plot_command,
										const char *build_name, int i,
										size_t *len,
										size_t plot_command_max_size,int nframes)
{
	
	int frames_bottleneck=MIN(nframes,series->ntimes/series->frameskip);
	
	int time = 0;
	
	if(series->is_static)
	{
		char data_path[100];
		snprintf(data_path,sizeof data_path, "data/%s_%d_%d_%05d.dat", 
			build_name, 
			gb_counter, i,0);
		double *carrs[2]={series->spatial_grid,series->values[0]};
		carrs_to_filename(data_path, carrs, series->npoints, 2);
	}
	else{
		for(int j=0; j<frames_bottleneck; j++)
		{
			char data_path[100];
			snprintf(data_path,sizeof data_path, "data/%s_%d_%d_%05d.dat", 
				build_name, 
				gb_counter, i,j);
			time = j*series->frameskip;
			double *carrs[2]={series->spatial_grid,series->values[time]};
			carrs_to_filename(data_path, carrs, series->npoints, 2);
		}
	}
	
	char *label=(series->label) ? series->label : "";
	char *style=(series->style) ? series->style : "lp";
	
	*len+=snprintf(plot_command+*len,plot_command_max_size-*len,
				"file%d using 1:2 with %s title \"%s\"",i,style, label);
}

void anim_set_terminal_and_plot(FILE *gp, char *plot_command, 
	const char *output_name,double delay)
{
	char script[5'000];
	char output_path[100];
	sprintf(output_path, "\'plot/%s.gif\'", output_name);
	sprintf(script,
		"set terminal gif animate delay %.2lf size 1000,700 font \"CMU Serif,16\" enhanced\n"
		"set for [i=1:8] linetype i lw 4.5 pt 7 ps 0.65\n"
		"set output %s\n"
		"%s\n"
		"unset output\n", delay,output_path, plot_command );
	fputs(script, gp);

}

void eplot_anim(const AnimSpec *series, int count, const Global *gb_settings,
			int nframes,const char *additional_in_loop,double delay)
{	
	char *build_name=get_build_name(gb_settings);
	
	char script_path[100]; 
	FILE *gp=set_and_open_script(build_name, script_path);
	
	eplot_set_common(gp);
	
	eplot_set_global(gp, gb_settings);
	
	char plot_command[2'000];
	size_t len=snprintf(plot_command,sizeof plot_command, "do for [i=0:%d] {\n",
		nframes-1);
		
	for(int i=0; i<count; i++)
	{
		anim_data_plt_command_helper1(plot_command,build_name,i,&len,
			sizeof plot_command,series[i].is_static);
	}
	
	len+=snprintf(plot_command+len,sizeof plot_command-len,"%s\nplot ",
		additional_in_loop);
	
	for(int i=0; i<count; i++)
	{	
		anim_data_plt_command_helper2(series+i, plot_command, build_name, i,&len,
			sizeof plot_command,nframes);
		len+=snprintf(plot_command+len,sizeof plot_command - len,
												(i==count-1)?"\n":", \\" "\n");
	}
	
	len+=snprintf(plot_command+len,sizeof plot_command - len,"}\n");
	
	eprint( "-----------------------------------------\n"
			"Creating animation '%s'\n"
			"Plot command written in .gp script:\n"
			"%s"
			"-----------------------------------------", 
			gb_settings->title,plot_command);
	anim_set_terminal_and_plot(gp, plot_command, gb_settings->output_name,delay);

	fflush(gp);
	fclose(gp);

	char shell_command[200];
	snprintf(shell_command, sizeof shell_command, "gnuplot %s", script_path);
	system(shell_command);
	
	puts("Plot was successful\n"
		 "-----------------------------------------\n");
	
	gb_counter++;
}

Global eplot_init_global(	
	char *title,
	char *xlabel,
	char *ylabel,
	char *logscale,
	char *build_name,
	char *output_name,
	char *key,
	char *additional_set)
{
	Global gb_settings={
		title,
		xlabel,
		ylabel,
		logscale,
		build_name,
		output_name,
		key,
		additional_set
	};
	
	return gb_settings;
}

SeriesSpec eplot_init_series(
	double *x,
	double *y,
	int size,
	char *style,
	char *label,
	char *color)
{
	SeriesSpec series ={
		x,
		y,
		size,
		style,
		label,
		color
	};
	return series;
}

AnimSpec eplot_init_animseries(
	double *spatial_grid,
	int npoints,
	double **values,
	int ntimes,
	int frameskip,
	char *style,
	char *label,
	char *color,
	bool is_static)
{
	AnimSpec series ={
		spatial_grid,
		npoints,
		values,
		ntimes,
		frameskip,
		style,
		label,
		color,
		is_static
	};
	return series;
}
