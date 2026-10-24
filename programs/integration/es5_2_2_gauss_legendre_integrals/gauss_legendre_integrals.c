#include "integration.h"
#include "explot.h"
#include <math.h>

#define MAX_N_NODES 40
#define STEP 4
#define MAX_LABEL_LENGTH 100

double f1( double x, void *unused ){
	UNUSED(unused);
	return exp(-4*x);}
	
double f2( double x, void *unused ){
	UNUSED(unused);
	return exp(-9*x*x);}
	
double f3( double x, void *unused ){
	UNUSED(unused);
	return 1/cosh(x);}

double f4( double x, void *unused ){
	UNUSED(unused);
	return 1/(1 + 9*x*x);}
	
double f5( double x, void *unused ){
	UNUSED(unused);
	return x*x*sin(8*x);}

void computation_result_analysis( int max_n_nodes, ParamFuncPtr f, double true_value,
	double x1, double x2, Array *true_value_deviation,char *fname)
{
	int n_nodes[ 1000 ];
	const int n_estimations = max_n_nodes / STEP;
	Array evaluated_integrals = arr_new(n_estimations);
	arr_init(true_value_deviation,n_estimations);
		
	printf( ANSI_GREEN"Evalutaing integral of function "ANSI_RED"%s\n"
			ANSI_GREEN"with Gauss-Legendre quadrature rule\n"
			ANSI_YELLOW"between endpoints [%.2lf,%.2lf]\n"
			"True value is "ANSI_RED"%.4lf\n\n",fname,x1,x2,true_value);
	
	puts(ANSI_YELLOW"n° nodes | diff with true value"ANSI_RESET);
	
	for(int i = 0; i<n_estimations; i++)
	{
		
		n_nodes[i] = ( i+1 ) * STEP;
		ARR(evaluated_integrals, i)=intg_gauss_legendre( f,NULL, n_nodes[i],
								x1, x2);
		ARRP(true_value_deviation, i)=fabs(true_value - ARR(evaluated_integrals,i));
		printf("%8d | %.3e\n",n_nodes[i],ARRP(true_value_deviation,i));
	}
	
	puts("\n-------------------------------------\n");

	arr_free(&evaluated_integrals);
}

void plot_error(const Array *true_value_deviation, int nfs,
			char fnames[][MAX_LABEL_LENGTH])
{
	const int n_estimations=MAX_N_NODES/STEP;
	Array n_nodes = arr_linspace(STEP,MAX_N_NODES,n_estimations);
	
	SeriesSpec *series;
	SAFE_ALLOC(series,nfs);
	
	for(int i=0; i<nfs; i++)
		series[i]=eplot_init_series(n_nodes.data,(true_value_deviation+i)->data,
				n_estimations,"lp",fnames[i],NULL);
	
	Global gb_settings = eplot_init_global("Convergence of Guass-Legendre"
		"quadrature rule","n. nodes","|I_n + I|","y","gauss_leg_conv_plot",
		"gauss_leg_conv_plot",NULL,
		"set object 1 rect from -20,1e-13 to 60,1e-16 fc rgb \"light-red\" " 
		"fs transparent solid 0.2 border lc rgb \"#b0f03232\" lw 1.5\n"
		"set object 1 front\n"
		"set label \"{/Symbol e}_{mach} zone\" at 1,4e-13\n"
			);
	
	eplot_multi_options(series, nfs, &gb_settings, false);
	
	arr_free(&n_nodes);
	SAFE_FREE(series);
}



int main()
{
	Array true_value_deviation[5]={0};
	char fnames[][MAX_LABEL_LENGTH]={
		{"e^{-4x}"},
		{"e^{-9x^2}"},
		{"sech(x)"},
		{"1/(1+9x^2)"},
		{"x^2sin(8x)"}
	};
	
	computation_result_analysis( MAX_N_NODES, f1, sinh(4)/2, -1, 1,
			true_value_deviation,fnames[0]);
	computation_result_analysis( MAX_N_NODES, f2, sqrt(PI) * erf(3)/3, -1, 1,
			true_value_deviation+1,fnames[1]);
	computation_result_analysis( MAX_N_NODES, f3, 2*atan( sinh(1) ), -1, 1,
			true_value_deviation+2,fnames[2]);
	computation_result_analysis( MAX_N_NODES, f4, 2*atan(3)/3, -1 ,1,
			true_value_deviation+3,fnames[3]);
	computation_result_analysis( MAX_N_NODES, f5, -3*PI*PI/32, PI/2 , PI,
			true_value_deviation+4,fnames[4]);
			
	plot_error(true_value_deviation,5,fnames);
			
	// Here it got the exact value in the second to last calculation
	// It's an artifact of floating point representation. Integrals
	// are all irrational and cannot be represented precisely on a computer
	
	for(int i=0; i<5; i++)
		arr_free(true_value_deviation+i);
	
	return 0;
}
