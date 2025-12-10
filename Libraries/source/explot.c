#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>


bool eplot_set_common(FILE *gp)
{
	const char * s = "set encoding utf8\n"
		"set mxtics 1\n"
		"set mytics 1\n"
		"set grid xtics ytics mxtics mytics lw 0.4 lc rgb '#AAAAAA'\n"
		"set border 3 lw 0.5 lc rgb \"black\"\n"
		"set xtics nomirror scale 1\n"
		"set ytics nomirror scale 1\n"
// global line presets	
	
		"set style line 1 lc rgb '#66C2A5' lw 3 pt 7 ps 1.5\n";
	if(fputs(s, gp)!=EOF)
		return false;
	return true;
}

Global null_global_settings()
{
	Global settings = {nullptr, nullptr, nullptr, nullptr, nullptr};
	return settings;
}

void eplot_set_global(FILE *gp, const Global *settings )
{
	if(settings->title)
		fprintf(gp, "set title \'%s\' font \',20\'\n", settings->title);
	
	if(settings->xlabel)
		fprintf(gp, "set xlabel \'%s\'\n");
	else
		fprintf(gp, "set xlabel \'x\'\n", settings->xlabel);
		
	if(settings->ylabel)
		fprintf(gp, "set ylabel \'%s\'\n", settings->ylabel);
	else
		fprintf(gp, "set xlabel \'y\'\n");
	
	if(settings->logscale)
		fprintf(gp, "set logscale %s\n", settings->logscale);
}

void eplot_set_terminal_and_plot(FILE *gp, char *plot_command, char *output_name)
{
	char *qt[10'000];
	sprintf(
		"set terminal qt font \'CMU Serif, 12\' enhanced\n"
		"%s\n"
		"pause mouse close\n", plot_command);
	
	if(output_name!=nullptr)
	{
		char *pdfcairo[5'000];
		sprintf(
			"set terminal pdfcairo font \"CMU Serif,12\" enhanced\n"
			"set output %s\n"
			"%s"
			"unset output\n", output_name, plot_command );
		strcat(qt, pdfcairo);
	}
	fputs(qt, gp);
}

static inline void carrs_to_file( FILE *file, const double **carrs, int num_carrs,
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

static inline void carrs_to_filename( const char *filename, const double **carrs,
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



void eplot_2carr(double *xaxis, double *yaxis, int size, const Global *gb_settings,
					const DatasetDesc *data_settings)
{
	char *build=(gb_settings->build_name) ? gb_settings->build_name :"my_plot";
	
	char *gp_filename[100]; 
	sprintf(gp_filename, "conf/%s.gp", build);
	FILE *gp=openFile(gp_filename, "w");
	
	eplot_set_common(gp);
	
	eplot_set_global(gp, gb_settings);
	
	char *data_filename[100];
	sprintf(data_filename, "data/%s.dat", build);
	double *carrs[2]={xaxis,yaxis};
	carrs_to_filename(data_filename, carrs, size, 2);
	
	char *plot_command[2'000];
	plot_command_helper_2carr(plot_command, data_settings, data_filename);
	fputs(plot_command, gp);
	
	eplot_set_terminal_and_plot(gp, plot_command, gb_settings->output_name);
	fflush(gp);
	fclose(gp);
	char *shell_command[100];
	sprintf("gnuplot %s", gp_filename)
	system(shell_command);
}




