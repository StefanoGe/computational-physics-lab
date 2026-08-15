#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rootf.h"
#include "explot.h"

#define DOMAIN_PLOT_SIZE 1000

#define ANSI_BLUE   "\x1b[34m"
#define ANSI_RESET  "\x1b[0m"
#define ANSI_RED    "\x1b[31m"
#define ANSI_YELLOW "\x1b[33m"
#define ANSI_GREEN "\x1b[32m"

typedef void(*c_value_printer)(ParamFuncPtr,FuncPtr,FuncPtr,double,const Array*);

void c_value_quadratic(ParamFuncPtr der, FuncPtr der2, FuncPtr unused,
						double root, const Array *abs_errs)
{
	Array c_ratios = arr_new(abs_errs->size-1);
	
	for(int i=0; i<c_ratios.size; i++)
		ARR(c_ratios,i) = ARRP(abs_errs,i+1) / pow(ARRP(abs_errs,i),2);
	
	UNUSED(unused);
	puts(ANSI_YELLOW"Ratios |abs_err(n+1)|/|abs_err(n)|^2\n"ANSI_RESET);
	arr_print_inline(&c_ratios,"%lf",true,true);
	putchar('\n');
	printf(ANSI_YELLOW"Expected C value: %lf\n\n"ANSI_RESET
				,fabs(der2(root)/der(root,NULL) )/2 );
	arr_free(&c_ratios);
}

void c_value_cubic(ParamFuncPtr der, FuncPtr unused, FuncPtr der3,
						double root, const Array *abs_errs)
{
	Array c_ratios = arr_new(abs_errs->size-1);
	
	for(int i=0; i<c_ratios.size; i++)
		ARR(c_ratios,i) = ARRP(abs_errs,i+1) / pow(ARRP(abs_errs,i),3);
	
	UNUSED(unused);
	puts(ANSI_YELLOW"Ratios |abs_err(n+1)|/|abs_err(n)|^2\n"ANSI_RESET);
	arr_print_inline(&c_ratios,"%lf",true,true);
	putchar('\n');
	printf(ANSI_YELLOW"Expected C value: %lf\n\n"ANSI_RESET,
					fabs(der3(root)/der(root,NULL) )/3 );
	arr_free(&c_ratios);
}


void c_value_linear(ParamFuncPtr unused1, FuncPtr unused2, FuncPtr unused3,
						double root, const Array *abs_errs)
{
	Array c_ratios = arr_new(abs_errs->size-1);
	
	for(int i=0; i<c_ratios.size; i++)
		ARR(c_ratios,i) = ARRP(abs_errs,i+1) / ARRP(abs_errs,i);
	
	UNUSED(unused1);
	UNUSED(unused2);
	UNUSED(unused3);
	UNUSED(root);
	puts(ANSI_YELLOW"Ratios |abs_err(n+1)|/|abs_err(n)|^2\n"ANSI_RESET);
	arr_print_inline(&c_ratios,"%lf",true,true);
	putchar('\n');
	printf(ANSI_YELLOW"Expected C value: %lf\n\n"ANSI_RESET,0.5 );
	arr_free(&c_ratios);
}

double f1( double x, void * unused )
{
	UNUSED(unused);
	return x*x - exp( -x );
}

double f1_der( double x, void * unused )
{
	UNUSED(unused);
	return 2 * x + exp(-x);
}

double f1_der2(double x){ return 2 - exp(-x); }

double f2( double x, void * unused )
{
	UNUSED(unused);
	return tan(x) - 2 * x;
}

double f2_der( double x, void * unused )
{
	UNUSED(unused);
	return 1/cos(x)/cos(x) - 2 ;
}

double f2_der2(double x){return -2*sin(x)/pow(cos(x),3);}

double f2_der3(double x){return -2*(pow(cos(x),2) 
							+ 3*pow(sin(x),2))/pow(cos(x),4);}

double f3( double x, void * unused )
{
	UNUSED(unused);
	return exp(x + 1) - x - 2;
}

double f3_der( double x, void * unused )
{
	UNUSED(unused);
	return exp(x+1) - 1;
}

double f3_der2(double x){return exp(x+1);}

void convergence_study(const Array *debug, const Array *roots_found, 
						int nfound_roots, const Array *starting_values,
						ParamFuncPtr der,FuncPtr der2,FuncPtr der3,
						const char *f_name,
						c_value_printer *c_printer)
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
		printf(ANSI_GREEN"Newton's method convergence study for root x = %lf\n"
				"of function "ANSI_RESET ANSI_RED"%s\n"ANSI_RESET
				ANSI_YELLOW"starting guess x = %.2lf\n\n"ANSI_RESET, 
			ARRP(roots_found,i),f_name,ARRP(starting_values,i));
		puts(ANSI_YELLOW"-log|abs_err|:"ANSI_RESET);
		arr_print_inline(log_abs_errs+i,"%lf",true,true);
		puts(ANSI_YELLOW"\nRatios of consecutive -log|abs_err|:"ANSI_RESET);
		arr_print_inline(log_ratios+i,"%.3lf",true,true);
		putchar('\n');
		(*c_printer[i])(der,der2,der3,ARRP(roots_found,i),abs_errs+i);
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

void test_newt(ParamFuncPtr f,ParamFuncPtr derivative,FuncPtr der2,
				FuncPtr der3,
				double x1,double x2,const double starting_points[], int nroots,
				char *f_name, c_value_printer *c_printer)
{
	Array *debug;
	SAFE_ALLOC(debug,nroots);
	for(int i=0; i<nroots; i++)
		arr_null(debug+i);
	Array roots_found = {0};
	Array starting_values={0};
	double curr_root = 0;
	int nfound_roots=0;
	
	for(int i=0; i<nroots; i++)
	{
		if( !isnan(curr_root=root_newt(f,NULL,derivative,NULL,starting_points[i],
								DEF_TOL,DEF_TOL,debug+nfound_roots)))
		{
			arr_append(&roots_found,curr_root);
			arr_append(&starting_values,starting_points[i]);
			nfound_roots++;
			for(int j=0; j<debug[nfound_roots-1].size; j++)
			{
				if(fabs(ARR(debug[nfound_roots-1],j) - 
				ARR(roots_found,nfound_roots-1) )< 1e-15 )
					debug[nfound_roots-1].size=j;
			}
		} else{
			arr_free(debug+nfound_roots);
			arr_null(debug+nfound_roots);
		}
	}
	
	convergence_study(debug,&roots_found,nfound_roots,&starting_values,
						derivative,der2,der3,f_name,c_printer);
	
	Array domain = arr_linspace(x1,x2,DOMAIN_PLOT_SIZE);
	Array y_values = arr_map_par(&domain,f,NULL);
	Array y_roots = arr_map_par(&roots_found,f,NULL);
	
	char title[100];
	snprintf(title,100,"Graph and roots of %s",f_name);
	
	Global gb_settings = eplot_init_global(title,
		"x","y",NULL,"test_newt",NULL,"top left",NULL);
	
	SeriesSpec graph = eplot_init_series(domain.data,y_values.data,
		domain.size,"l",f_name,NULL);
	
	SeriesSpec roots = eplot_init_series(roots_found.data,y_roots.data,
		roots_found.size,"p","Roots",NULL);
	
	SeriesSpec series[]={graph,roots};
	
	eplot_multi(series,2,&gb_settings);
	putchar('\n');
	
	for(int i=0; i<nfound_roots; i++)
		arr_free(debug+i);
	SAFE_FREE(debug);
	Array *arr_to_free[]={&domain,&y_values,&roots_found,&y_roots,&starting_values};
	arr_free_many(arr_to_free,5);
}

int main()
{
	const double start1[]={-1};
	c_value_printer print_funcs1[]={c_value_quadratic};
	test_newt(f1,f1_der,f1_der2,NULL,-2,2,start1,1,"x^2 - e^{-x}",
		print_funcs1);
//  We find quadratic convergence. Asymptotic constant C is what is expected.
	
	c_value_printer print_funcs2[]={c_value_cubic,c_value_quadratic};
	const double start2[]={-0.2,1,1.4};
	test_newt(f2,f2_der,f2_der2,f2_der3,-0.2,1.4,start2,3,"tan(x) - 2x",
				print_funcs2);
//	2nd has quadratic convergence as expected. First root has second derivative
//	equal to 0 and third derivative different than 0 --> cubic convergence!
	
	c_value_printer print_funcs3[]={c_value_linear,c_value_linear,c_value_linear,
						c_value_linear};
	const double start3[]={-2,0,1,2};
	test_newt(f3,f3_der,NULL,NULL,-2,2,start3,4,"e^{x+1} - x - 2",print_funcs3);
//  Here convergence is slower. We are able to prove that it is linear with
//  asymptotic constant equal to 0.5.

	exit(EXIT_SUCCESS);
}


