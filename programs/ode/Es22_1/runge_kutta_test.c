#include "ode.h"
#include "explot.h"
#include "array.h"
#include "matrix.h"
#include <math.h>
#include "genutil.h"

double eq1(double t, double *u, [[maybe_unused]]void *unused ){
	return -2*t**u;}

double true1(double t){
	return 2*exp(-t*t);}

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
//	mypause();
	return max;
}

void error_study(Derivative eq, double t0, double tf, double in_value,
	const char *fnc_name, FuncPtr fnc)
{
	const int max_nstep=300;
	const int nstep_increment=30;
	const int n_iter=max_nstep/nstep_increment;
	// Variables passed to ode euler solver
	Matrix num_sol={0};
	Array time={0};
	// Global errors
	Array errors_ie2={0};
	Array errors_rk4={0};
	// Number of function evaluation for each iteration of the methods
	Array nf_ie2={0};
	Array nf_rk4={0};

	Array true_values={0};
	Array nsteps_vec={0};
	Array one_over_step_size_sq={0};
	Array one_over_step_size_4th={0};
	
	for(int i=0; i<n_iter; i++)
	{
		// Number of steps
		int nsteps=nstep_increment*(i+1);
		arr_append(&nf_ie2, (double)nsteps*2);
		arr_append(&nf_rk4, (double)nsteps*4);
		arr_append(&one_over_step_size_sq, (double)1/nsteps/nsteps);
		arr_append(&one_over_step_size_4th, ARR(one_over_step_size_sq,i)/nsteps/nsteps);
		arr_append(&nsteps_vec, (double)nsteps);
		// Call of IE2 method
		ode_ie2( t0, tf, &in_value, &eq, NULL,
						1, nsteps, &time, &num_sol );
		// Computation of true values
		arr_init(&true_values, nsteps+1);
		for(int j=0; j<nsteps+1; j++)
			ARR(true_values,j)=fnc(ARR(time,j));
		
		// Error computation for IE2 method
		arr_append(&errors_ie2,  compute_error(true_values.data, num_sol.data,nsteps+1) );
		
		// Call of RK4 method
		ode_rk4( t0, tf, &in_value, &eq, nullptr,
						1, nsteps, &time, &num_sol );
		
		//Error computation for RK4 method
		arr_append(&errors_rk4, compute_error(true_values.data, num_sol.data,nsteps+1) );		
		arr_free(&time);
		mat_free(&num_sol);
	}
	// Plot setup --- global settings
	char title[1000];
	sprintf(title, "Convergence of IE2 and RK4 method for u\'(t)=%s", 
							fnc_name);

	Global gb_settings = eplot_init_global(title,"Number of function evaluations",
		"Maximum error","xy","ie2_rk4_conv_study","ie2_rk4_conv_study",NULL,NULL);
	// Max error over whole vs nsteps plot 
	SeriesSpec error_ie2=eplot_init_series(nf_ie2.data,errors_ie2.data,n_iter,
		"lp","ie2",NULL);
	// Max error over whole vs nsteps plot 
	SeriesSpec error_rk4=eplot_init_series(nf_rk4.data,errors_rk4.data,
		n_iter,"lp","rk4",NULL);
	// Reference line for second order convergence
	SeriesSpec ref_line2={nsteps_vec.data,one_over_step_size_sq.data,n_iter,
		"l dt 4","2^{nd} order convergence",NULL};
	// Reference line for fourth order convergence
	SeriesSpec ref_line4={nsteps_vec.data,one_over_step_size_4th.data,
		n_iter,"l dt 4","4^{th} order convergence",NULL};
	SeriesSpec series[]={error_ie2, error_rk4, ref_line2, ref_line4};
	eplot_multi(series, 4, &gb_settings);
			
	Array *arr_to_free[]={&errors_ie2,&errors_rk4,&nf_ie2,&nf_rk4,
		&true_values,&nsteps_vec,&one_over_step_size_4th,&one_over_step_size_sq};
	arr_free_many(arr_to_free,8);
}

int main()
{
	error_study(eq1, 0, 2, 2, "-2tu", true1);
	
	return 0;
}
