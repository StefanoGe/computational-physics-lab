#include "integration.h"
#include <math.h>
#include <limits.h>
#include <float.h>
#include "explot.h"

#define MAXN 60
#define MINN 4
#define STEPN 2
#define MAX_LABEL_LENGTH 100
#define NTESTS 5

const int n_estimations = (MAXN-MINN) / STEPN +1;

double f1(double x, void *unused){
			UNUSED(unused);
	return 1/(1+x*x+x*x*x*x);}

double f2(double x, void *unused){
	UNUSED(unused);
	return exp(-x*x)*cos(x);}
	
double f3(double x, void *unused){
	UNUSED(unused);
	return pow(1.0+x*x,-2.0/3.0);}
	
double f4(double x, void *unused){
	UNUSED(unused);
	return 1.0/(1+x*x);}

double f5(double x, void *unused){
	UNUSED(unused);
	return exp(-x)/sqrt(x);}

void computation_result_analysis( ParamFuncPtr f, double true_value,
			Array *true_value_deviation,char *fname, Dbexp_mode mode,
			char *endpoint1,char *endpoint2)
{
	int n_nodes[ 1000 ];
	Array evaluated_integrals = arr_new(n_estimations);
	arr_init(true_value_deviation,n_estimations);
	
	printf( ANSI_GREEN"Evalutaing integral of function "ANSI_RED"%s\n"
			ANSI_GREEN"with double exponential quadrature rule\n"
			ANSI_YELLOW
			"between endpoints ["ANSI_RED"%s"ANSI_RESET","ANSI_RED"%s"ANSI_RESET
			"]\n"
			"True value is "ANSI_RED"%.4lf\n\n",fname,
			endpoint1,endpoint2,true_value);
	
	puts(ANSI_YELLOW"n° nodes | diff with true value"ANSI_RESET);
	for(int i=0; i<n_estimations; i++)
	{
		n_nodes[i] = MINN+i*STEPN;
		ARR(evaluated_integrals,i)=intg_double_exp(f, NULL, n_nodes[i]/2, mode,
				0,0, DEF_TOL);
		ARRP(true_value_deviation, i)=fabs(true_value - ARR(evaluated_integrals,i));
		printf("%8d | %+.3e\n",n_nodes[i], true_value - ARR(evaluated_integrals,i));
	}
	puts("\n-------------------------------------\n");
	arr_free(&evaluated_integrals);
}

void plot_error(const Array *true_value_deviation, int nfs,
			char fnames[][MAX_LABEL_LENGTH])
{
	Array n_nodes = arr_linspace(MINN,MAXN,n_estimations);
	
	SeriesSpec *series;
	SAFE_ALLOC(series,nfs);
	
	for(int i=0; i<nfs; i++)
		series[i]=eplot_init_series(n_nodes.data,(true_value_deviation+i)->data,
				n_estimations,"lp",fnames[i],NULL);
	
	Global gb_settings = eplot_init_global("Convergence of double exp. "
		"rule for improper integrals",
		"n. nodes","|I_n - I|","y","improper_conv_plot",
		"improper_conv_plot",NULL,
		"set object 1 rect from -20,1e-14 to 60,1e-16 fc rgb \"light-red\" " 
		"fs transparent solid 0.2 border lc rgb \"#b0f03232\" lw 1.5\n"
		"set object 1 front\n"
		"set label \"{/Symbol e}_{mach} zone\" at 1,4e-14\n"
			);
	
	eplot_multi_options(series, nfs, &gb_settings, false);
	
	arr_free(&n_nodes);
	SAFE_FREE(series);
}

int main()
{
	Array true_value_deviation[NTESTS]={0};
	char fnames[][MAX_LABEL_LENGTH]={
		{"1/(1+x^2+x^4)"},
		{"e^{-x^2}cos(x)"},
		{"(1+x^2)^{-2/3}"},
		{"1/(1+9x^2)"},
		{"e^{-x}/sqrt(x)"}
	};
	
	computation_result_analysis(f1,PI/sqrt(3),true_value_deviation,fnames[0],
								REAL_LINE,"-inf","+inf");
	computation_result_analysis(f2,sqrt(PI)/exp(0.25),
								true_value_deviation+1,fnames[1],REAL_LINE
								,"-inf","+inf");
	computation_result_analysis(f3,sqrt(PI)*tgamma(1.0/6)/tgamma(2.0/3),
								true_value_deviation+2,fnames[2],REAL_LINE
								,"-inf","+inf");
	computation_result_analysis(f4,PI/2,true_value_deviation+3,
								fnames[3],HALF_LINE
								,"0","+inf");
	computation_result_analysis(f5,sqrt(PI),true_value_deviation+4,
								fnames[4],HALF_LINE_EXP
								,"0","+inf");
	
	plot_error(true_value_deviation,5,fnames);
	
	for(int i=0; i<NTESTS; i++)
		arr_free(true_value_deviation+i);
	
	exit(EXIT_SUCCESS);
}
