
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
	return sqrt(x)*log(x);}

double f2(double x, void *unused){
	UNUSED(unused);
	return sqrt(1-x*x);}
	
double f3(double x, void *unused){
	UNUSED(unused);
	if(x==0)
		raiseErr("cannot take log(0)");
	return pow(log(x),2);}
	
double f4(double x, void *unused){
	UNUSED(unused);
	return log(cos(x));}

double f5(double x, void *unused){
	UNUSED(unused);
	if(x==PI/2)
		return NAN;
	return sqrt(tan(x));}

void computation_result_analysis( ParamFuncPtr f, double true_value,
			Array *true_value_deviation_dbe,Array *true_value_deviation_gl,
			char *fname,
			double endpoint1,double endpoint2,double tol)
{
	int n_nodes[ 1000 ];
	Array evaluated_integrals_dbe = arr_new(n_estimations);
	Array evaluated_integrals_gl = arr_new(n_estimations);
	arr_init(true_value_deviation_dbe,n_estimations);
	arr_init(true_value_deviation_gl,n_estimations);
	
	printf( ANSI_GREEN"Evalutaing integral of function "ANSI_RED"%s\n"
			ANSI_GREEN"with double exponential quadrature rule\n"
			ANSI_YELLOW
			"between endpoints ["ANSI_RED"%.2lf"ANSI_RESET","ANSI_RED
			"%.2lf"ANSI_RESET
			"]\n"
			"True value is "ANSI_RED"%.4lf\n\n",fname,
			endpoint1,endpoint2,true_value);
	
	puts(ANSI_YELLOW"n° nodes | diff with true value");
	for(int i=0; i<11;i++)
		putchar(' ');
	puts("db exp    | gauss-legendre"ANSI_RESET);
	for(int i=0; i<n_estimations; i++)
	{
		n_nodes[i] = MINN+i*STEPN;
		ARR(evaluated_integrals_dbe,i)=intg_double_exp( f,NULL, n_nodes[i],BOUNDED,
			endpoint1,endpoint2,tol);
		ARR(evaluated_integrals_gl,i)=intg_gauss_legendre(f,NULL,n_nodes[i],
			endpoint1,endpoint2);
		ARRP(true_value_deviation_dbe, i)=fabs(true_value - ARR(evaluated_integrals_dbe,i));
		ARRP(true_value_deviation_gl,i)=fabs(true_value - ARR(evaluated_integrals_gl,i));
		printf("%8d | %.3e | %.3e\n",n_nodes[i],ARRP(true_value_deviation_dbe,i),
				ARRP(true_value_deviation_gl,i));
	}
	puts("\n-------------------------------------\n");
	arr_free(&evaluated_integrals_dbe);
	arr_free(&evaluated_integrals_gl);
}

void plot_error(const Array *true_value_deviation_dbe, 
			const Array *true_value_deviation_gl, int nfs,
			char fnames[][MAX_LABEL_LENGTH])
{
	Array n_nodes = arr_linspace(MINN,MAXN,n_estimations);
	
	SeriesSpec *series;
	SAFE_ALLOC(series,nfs*2);
	
	char fnames_de[NTESTS][MAX_LABEL_LENGTH];
	char linetype_gl[NTESTS][MAX_LABEL_LENGTH];
	char linetype_de[NTESTS][MAX_LABEL_LENGTH];
	
	for(int i=0; i<nfs; i++)
	{
		sprintf(fnames_de[i],"%s de", fnames[i]);
		sprintf(linetype_gl[i],"lp lt %d",i+1);
		sprintf(linetype_de[i],"l lt %d dt 2",i+1);
		series[i*2]=eplot_init_series(n_nodes.data,(true_value_deviation_dbe+i)->data,
				n_estimations,linetype_de[i],fnames_de[i],NULL);
		series[i*2+1]=eplot_init_series(n_nodes.data,(true_value_deviation_gl+i)->data,
				n_estimations,linetype_gl[i],"-    gl",NULL);
	}

	
	Global gb_settings = eplot_init_global("Convergence of double exp. "
		"rule for improper integrals",
		"n. nodes","|I_n - I|","y","improper_conv_plot",
		"improper_conv_plot","outside",
		"set object 1 rect from -20,1e-14 to 60,1e-16 fc rgb \"light-red\" " 
		"fs transparent solid 0.2 border lc rgb \"#b0f03232\" lw 1.5\n"
		"set object 1 front\n"
		"set label \"{/Symbol e}_{mach} zone\" at 1,4e-14\n"
//		"set xrange [0:80]\n"
			);
	
	eplot_multi(series,nfs*2,&gb_settings);
	
	arr_free(&n_nodes);
	SAFE_FREE(series);
}

int main()
{
	// Problem: we need to tabulate some values to understand and decide
	// how to choose tM
/*
	printf("Maximum double: %e\n",DBL_MAX);

	for(int i=0; i<100; i++)
	{
		printf("%.2lf %e\n",i*0.1,sinh(PI/2*sinh(i*0.1)));
	}
*/
	Array true_value_deviation_dbe[NTESTS]={0};
	Array true_value_deviation_gl[NTESTS]={0};
	char fnames[][MAX_LABEL_LENGTH]={
		{"sqrt(x)log(x)"},
		{"sqrt(1-x^2)"},
		{"(log(x))^2"},
		{"log(cos(x))"},
		{"sqrt(tan(x))"}
	};

	computation_result_analysis(f1,-4.0/9,true_value_deviation_dbe,
								true_value_deviation_gl,fnames[0],
								0,1,DEF_TOL);
	computation_result_analysis(f2,PI/4,
								true_value_deviation_dbe+1,
								true_value_deviation_gl+1,fnames[1],0,1,
								DEF_TOL);
	computation_result_analysis(f3,2,
								true_value_deviation_dbe+2,
								true_value_deviation_gl+2,fnames[2],0,1,
								10e-10);
	computation_result_analysis(f4,-PI/2*log(2),true_value_deviation_dbe+3,
								true_value_deviation_gl+3,
								fnames[3],0,PI/2,DEF_TOL);
	computation_result_analysis(f5,PI/sqrt(2),true_value_deviation_dbe+4,
								true_value_deviation_gl+4,
								fnames[4],0,PI/2,DEF_TOL);

//	printf("%e\n",(PI/2-10e-10)-(PI/2));

	plot_error(true_value_deviation_dbe,true_value_deviation_gl,5,fnames);
	
	for(int i=0; i<NTESTS; i++){
		arr_free(true_value_deviation_dbe+i);
		arr_free(true_value_deviation_gl+i);
	}
	
	exit(EXIT_SUCCESS);
}
