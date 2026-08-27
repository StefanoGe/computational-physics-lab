#include "integration.h"
#include <math.h>
#include <limits.h>
#include <float.h>

#define MAXN 60
#define MINN 4
#define STEPN 2
#define MAX_LABEL_LENGTH 100

double f1(double x, void *unused){
			UNUSED(unused);
	return 1/(1+x*x+x*x*x*x);}

double f2(double x, void *unused){
	UNUSED(unused);
	return exp(-x*x)*cos(x);}
	
double f3(double x, void *unused){
	UNUSED(unused);
	return pow(1.0+x*x,-2.0/3.0);}

void computation_result_analysis( ParamFuncPtr f, double true_value,
			Array *true_value_deviation,char *fname)
{
	int n_nodes[ 1000 ];
	const int n_estimations = (MAXN-MINN) / STEPN +1;
	Array evaluated_integrals = arr_new(n_estimations);
	arr_init(true_value_deviation,n_estimations);
	
	printf( ANSI_GREEN"Evalutaing integral of function "ANSI_RED"%s\n"
			ANSI_GREEN"with Gauss-Legendre quadrature rule\n"
			ANSI_YELLOW
			"True value is "ANSI_RED"%.4lf\n\n",fname,true_value);
	
	puts(ANSI_YELLOW"n° nodes | diff with true value"ANSI_RESET);
	for(int i=0; i<n_estimations; i++)
	{
		n_nodes[i] = MINN+i*STEPN;
		ARR(evaluated_integrals,i)=intg_double_exp_real_line( f,NULL, n_nodes[i]);
		ARRP(true_value_deviation, i)=fabs(true_value - ARR(evaluated_integrals,i));
		printf("%8d | %.3e\n",n_nodes[i],ARRP(true_value_deviation,i));
	}
	

	
	puts("\n-------------------------------------\n");


	arr_free(&evaluated_integrals);
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
	Array true_value_deviation[5]={0};
	char fnames[][MAX_LABEL_LENGTH]={
		{"1/(1+x^2+x^4)"},
		{"e^{-x^2}cos(x)"},
		{"(1+x^2)^{-2/3}"}
	};
	
	computation_result_analysis(f1,PI/sqrt(3),true_value_deviation,fnames[0]);
	computation_result_analysis(f2,sqrt(PI)/exp(0.25),
								true_value_deviation+1,fnames[1]);
	computation_result_analysis(f3,sqrt(PI)*tgamma(1.0/6)/tgamma(2.0/3),
								true_value_deviation+2,fnames[2]);

	
	for(int i=0; i<5; i++)
		arr_free(true_value_deviation+i);
	
	exit(EXIT_SUCCESS);
}
