#include "comp_physics.h"
#include "ode.h"
#include "explot.h"
#include <math.h>

void mypause ( void ) 
{ 
  printf ( "Press [Enter] to continue . . ." );
  fflush ( stdout );
  getchar();
} 

double eq1(double t, double *u, [[maybe_unused]]void *unused )
{return -2*t**u;}

double eq2(double t, double *u, [[maybe_unused]]void *unused )
{return *u+t;}

double eq3(double t, double *u, [[maybe_unused]]void *unused )
{return t*t/(*u)/(1+t*t*t);}

int powint(int base, int exp)
{
	int prod=1;
	for(int i=0; i<exp; i++)
		prod*=base;
	return prod;
}

typedef double(*TrueFunc)(double);

double true1(double t)
{
	return 2*exp(-t*t);
}

double true2(double t)
{
	return -1-t+3*exp(t);
}

double true3(double t)
{
	return sqrt(1 + 2*log(1+t*t*t)/3);
}

void solve_and_plot(Derivative eq, double t0, double tf, double in_value, 
		int nsteps, const char *fnc_name, TrueFunc func)
{
	Matrix f_output={0};
	Array t_output={0};
	Array initial_value=arr_asarr(&in_value, 1);
	static int counter=0;
	counter++;	

	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "Numerical solution of du/dt(t) = %s, u(%.0f)=%.0lf", 
							fnc_name, t0, in_value);
	sprintf(build_name, "euler_test%d", counter);
	sprintf(output_name, "euler_test%d", counter);

	Global gb_settings = {
		.title=title,
		.xlabel="t",
		.ylabel="u(t)",
		.logscale=nullptr,
		.build_name=build_name,
		.output_name=output_name
	};

	ode_feuler( t0, tf, &initial_value, &eq, 1, nsteps, &t_output, &f_output );

	Array true_values={0};
	arr_init(&true_values, nsteps+1);
	for(int i=0; i<nsteps+1; i++)
		true_values.data[i]=func(t_output.data[i]);

	SeriesSpec data_desc_true={
		.x=t_output.data,
		.y=true_values.data,
		.size=nsteps+1,
		.style="lines",
		.label="Exact function",
		.color=nullptr
	};
	
	SeriesSpec data_desc_est={
		.x=t_output.data,
		.y=f_output.data,
		.size=nsteps+1,
		.style="points",
		.label="Numerical solution",
		.color=nullptr
	};
	
	SeriesSpec specs[2]={data_desc_est, data_desc_true};

	eplot_multi(specs, 2, &gb_settings);


	
	mat_free(&f_output);
	arr_free(&t_output);
	arr_free(&true_values);
}

double compute_error(double *true_arr, double *estimate, int length)
{
	double max=0;
	double temp=0;
	for(int i=0; i<length; i++)
	{
		temp=fabs(true_arr[i]-estimate[i]);
		if(temp>=max)
			max=temp;
//		eprint("True value: %.6lf ; Estimate: %.6lf; max = %lf", true_arr[i], estimate[i], max);
	}
	return max;
}

void error_study(Derivative eq, double t0, double tf, double in_value,
	const char *fnc_name, int max_k, TrueFunc fnc)
{
	Matrix f_output={0};
	Array t_output={0};
	const Array initial_value=arr_asarr(&in_value, 1);
	VectorD errors=init_vecD();
	VectorD nsteps_vec=init_vecD();
	static int counter=0;
	counter++;
	Array true_values={0};
	
	for(int i=0; i<=max_k-2; i++)
	{
		int nsteps=10*powint(2,i+2);
		appendD(&nsteps_vec, (double)nsteps);
		
		ode_feuler( t0, tf, &initial_value, &eq, 1, nsteps, &t_output, &f_output );
		
		// Error computation
		arr_init(&true_values, nsteps+1);
		for(int j=0; j<nsteps+1; j++)
			true_values.data[j]=fnc(t_output.data[j]);
		appendD(&errors,compute_error(true_values.data, f_output.data,nsteps+1) );
		
		arr_free(&t_output);
		mat_free(&f_output);
	}
	
	// Plot setup
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "Convergence of Euler method for u(t)=%s, %.0lf=%.0lf", 
							fnc_name, t0, in_value);
	sprintf(build_name, "euler_conv_study%d", counter);
	sprintf(output_name, "euler_conv_study%d", counter);

	Global gb_settings = {
		.title=title,
		.xlabel="Number of steps",
		.ylabel="Absolute error",
		.logscale="xy",
		.build_name=build_name,
		.output_name=output_name
	};

	DatasetDesc data_desc={
		.file=nullptr,
		.style="lp",
		.label="",
		.color=nullptr
	};

	eplot_2carr(nsteps_vec.val, errors.val, max_k-1, &gb_settings,
			&data_desc);
			
	free_vecD(&errors, &nsteps_vec);
}

	
void test_feuler(Derivative eq, double t0, double tf, double in_value, 
	const char *fnc_name, TrueFunc true_func)
{
	solve_and_plot(eq, t0, tf, in_value, 320, fnc_name, true_func);
	error_study(eq, t0, tf, in_value, fnc_name, 10, true_func);
}


int main()
{
	
	test_feuler(eq1, 0, 2, 2, "-2tu", true1);
	test_feuler(eq2, 0, 1, 2, "u+t", true2);
	test_feuler(eq3, 0, 3, 1, "", true3);
	
	return 0;
}

/* OLD
	[[maybe_unused]]DatasetDesc data_desc={
		.file=nullptr,
		.style="lines",
		.label="",
		.color=nullptr
	};
//	eplot_2carr(t_output.data, f_output.data, nsteps+1, &gb_settings,
//			&data_desc);

*/
