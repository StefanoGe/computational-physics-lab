// Secant method implementation and study of convergence

#include <stdio.h>
#include <stdlib.h>
#include "rootf.h"
#include <math.h>
#include "explot.h"

#define DOMAIN_PLOT_SIZE 1000

#define ANSI_BLUE   "\x1b[34m"
#define ANSI_RESET  "\x1b[0m"
#define ANSI_RED    "\x1b[31m"
#define ANSI_YELLOW "\x1b[33m"
#define ANSI_GREEN "\x1b[32m"

double f1( double x, void * unused )
{
	UNUSED(unused);
	return x*x - exp( -x );
}

double f2( double x, void * unused )
{
	UNUSED(unused);
	return tan(x) - 2 * x;
}

double f3( double x, void * unused )
{
	UNUSED(unused);
	return exp(x + 1) - x - 2;
}


void plot_graph(ParamFuncPtr f, double x1, double x2, const Array *roots,
						char *fname)
{
	Array domain = arr_linspace(x1,x2,DOMAIN_PLOT_SIZE);
	Array y_values = arr_map_par(&domain,f,NULL);
	Array y_roots = arr_map_par(roots,f,NULL);
	
	char title[100];
	snprintf(title,100,"Graph and roots of %s",fname);
	
	Global gb_settings = eplot_init_global(title,
		"x","y",NULL,"test_newt",NULL,"top left",NULL);
	
	SeriesSpec graph = eplot_init_series(domain.data,y_values.data,
		domain.size,"l",fname,NULL);	

	SeriesSpec root_graph = eplot_init_series(roots->data,y_roots.data,
		roots->size,"p","Roots",NULL);

	SeriesSpec series[]={graph,root_graph};
	
	eplot_multi_options(series, 2, &gb_settings, false);
}

void convergence_study(const Array *debug, const Array *roots_found, 
						int nfound_roots, double *starting_values1,
						double *starting_values2,
						const char *f_name)
{
	Array *log_abs_errs;
	SAFE_ALLOC(log_abs_errs,nfound_roots);
	Array *log_ratios;
	SAFE_ALLOC(log_ratios,nfound_roots);
	Array *abs_errs;
	SAFE_ALLOC(abs_errs,nfound_roots);
	
	for(int root_index=0; root_index<nfound_roots; root_index++)
	{
		log_abs_errs[root_index] = arr_new(debug[root_index].size);
		abs_errs[root_index] = arr_new(debug[root_index].size);

		log_ratios[root_index] = arr_new(debug[root_index].size-1);
		for(int err_index=0; err_index<debug[root_index].size; err_index++)
		{
			ARR(abs_errs[root_index],err_index) = 
				fabs(ARRP(debug+root_index,err_index) 
				- ARRP(roots_found,root_index)) ;
			
			ARR(log_abs_errs[root_index],err_index) = -log10(
				ARR(abs_errs[root_index],err_index));
		}
			
		for(int err_index=0; err_index<debug[root_index].size-1; err_index++)
			ARR(log_ratios[root_index],err_index) = 
				ARR(log_abs_errs[root_index],err_index+1) /
				ARR(log_abs_errs[root_index],err_index);
	}
	
	puts("-----------------------------------\n");
	for(int i=0; i<nfound_roots; i++)
	{
		printf(ANSI_GREEN"Secant method convergence study for root x = %lf\n"
				"of function "ANSI_RESET ANSI_RED"%s\n"ANSI_RESET
				ANSI_YELLOW"starting interval = [%.2lf,%.2lf]\n\n"ANSI_RESET, 
			ARRP(roots_found,i),f_name,starting_values1[i],starting_values2[i]);
		puts(ANSI_YELLOW"-log|abs_err|:"ANSI_RESET);
		arr_print_inline(log_abs_errs+i,"%lf",true,true);
		puts(ANSI_YELLOW"\nRatios of consecutive -log|abs_err|:"ANSI_RESET);
		arr_print_inline(log_ratios+i,"%.3lf",true,true);
		putchar('\n');
		// (*c_printer[i])(der,der2,der3,ARRP(roots_found,i),abs_errs+i);
		if(i!=nfound_roots-1)
			puts("---\n");
	}
	
	puts("-----------------------------------\n\n");
	
	for(int i=0; i<nfound_roots; i++)
	{
		arr_free(log_ratios+i);
		arr_free(log_abs_errs+i);
		arr_free(abs_errs+i);
	}
	SAFE_FREE(log_ratios);
	SAFE_FREE(log_abs_errs);
	SAFE_FREE(abs_errs);
}

void remove_roots_from_debug(Array *debug, int nfound_roots,
		const Array *roots_found)
{
	for(int j=0; j<debug[nfound_roots-1].size; j++)
	{
		if(fabs(ARR(debug[nfound_roots-1],j) - 
		ARRP(roots_found,nfound_roots-1) )< 1e-15 )
			debug[nfound_roots-1].size=j;
	}
}

void test( ParamFuncPtr f, double x1, double x2,
					double xtol, double ftol, int nroots, 
					double starting_points1[], double starting_points2[],
					char *fname)
{
	Array *debug;
	SAFE_ALLOC(debug,nroots);
	for(int i=0; i<nroots; i++)
		arr_null(debug+i);
	Array roots = arr_new(nroots);
	
	for( int i = 0; i < nroots; i++)
	{
		int exit_status = -2;
		ARR(roots,i) = root_sec( f, NULL, starting_points1[i], 
							starting_points2[i], xtol, ftol, debug+i,
							&exit_status);
		
		if(exit_status == -1)
		{
			eprint("root_finding failed. debug:\n");
			arr_print( stderr,debug+i, "%.20lf, ");
			return;
		}
		remove_roots_from_debug(debug,i+1,&roots);
	}
	 
	convergence_study(debug,&roots,nroots,starting_points1,starting_points2,
		fname);
	 
	plot_graph( f, x1, x2, &roots,fname );

	arr_free(&roots);
	for(int i=0; i<nroots; i++)
		arr_free(debug+i);
	SAFE_FREE(debug);
}

int main()
{
	double starting_points1_1[]={-2};
	double starting_points1_2[]={2};
	test(f1,-2,2,DEF_TOL,DEF_TOL,1,starting_points1_1,starting_points1_2,
		"x^2 - e^{-x}");
	
	double starting_points2_1[]={-0.2,0.9};
	double starting_points2_2[]={0.5,1.4};
	test(f2,-0.2,1.4,DEF_TOL,DEF_TOL,2,starting_points2_1,starting_points2_2,
		"2x - tan(x)");
		
	double starting_points3_1[]={-2};
	double starting_points3_2[]={2};
	test(f3,-2,2,DEF_TOL,DEF_TOL,1,starting_points3_1,starting_points3_2,
		"e^{x+1} - 2 - x");

	exit(EXIT_SUCCESS);
}

