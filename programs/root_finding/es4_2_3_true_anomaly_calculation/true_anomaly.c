#include "rootf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "explot.h"

#define EROS_TAU 1.761
#define EROS_ECC 0.223
#define N_VALUES_TO_PLOT 100

const double two_pi_over_tau = 2*PI/EROS_TAU;

double positive_fmod(double x, double y)
{
	if(y<=0)
		raiseErr("y must be greater than 0");

	int int_division = (int) (x/y);
	
	if(x<0)
	{
		return x - (int_division - 1)*y;
	}
	else
	{
		return x - int_division*y;
	}
}

/*
 * params[0] = time
 * params[1] = eccentricity
 * params[2] = orbital period
 * 
 */

double ecc_an_eq(double psi, void *params)
{
	return psi - ((double*)params)[1]*sin(psi) 
		- 2*PI*((double*)params)[0]/((double*)params)[2];
}

double ecc_an_eq_der(double psi, void *params)
{
	return 1 - ((double*)params)[1] * cos(psi);
}

double true_anomaly(double psi, double eps)
{
	return 2*atan( sqrt((1+eps)/(1-eps)) * tan(psi/2) );
}

void true_anomaly_over_time()
{
	Array times = arr_linspace(0,EROS_TAU,N_VALUES_TO_PLOT);
	Array true_anomalies = arr_new(N_VALUES_TO_PLOT);
	
	double params[]= {0,EROS_ECC,EROS_TAU};
	
	for(int i=0; i<N_VALUES_TO_PLOT; i++)
	{
		params[0]=ARR(times,i);
		const double init_guess2 = two_pi_over_tau*ARR(times,i) - EROS_ECC;
		const double init_guess1 = two_pi_over_tau*ARR(times,i) + EROS_ECC;
		ARR(true_anomalies,i)=root_newt_bracket(ecc_an_eq, params, 
			ecc_an_eq_der, params, init_guess1, init_guess2,
			DEF_TOL, DEF_TOL, NULL);
		ARR(true_anomalies,i) = positive_fmod(true_anomaly(
			ARR(true_anomalies,i),EROS_ECC), 2*PI);
	}
	
	Global gb_settings = eplot_init_global("Plot of true anomaly {/Symbol q}(t)"
						" of asteroid Eros",
						"t","{/Symbol q}(t)",NULL,"asteroid_eros_plot",NULL,
						"top left",NULL);
	
	SeriesSpec plot_series = eplot_init_series(times.data,true_anomalies.data,
		times.size,"l","{/Symbol q}(t)",NULL);
	
	
	eplot_multi_options(&plot_series,1,&gb_settings, false);
	
	arr_free(&times);
	arr_free(&true_anomalies);
}


void simple_graph(double x1, double x2)
{
	double params[]= {0,EROS_ECC,EROS_TAU};
	Array times = arr_linspace(x1,x2,N_VALUES_TO_PLOT);
	Array y_values_eq = arr_map_par(&times,ecc_an_eq,params);
	
	Global gb_settings = eplot_init_global("","","",NULL,"equation_graph",
		NULL,"",NULL);
	
	SeriesSpec series = eplot_init_series(times.data,y_values_eq.data,times.size,
			"l",NULL,NULL);

	eplot_multi(&series,1,&gb_settings);
	
	arr_free(&times);
	arr_free(&y_values_eq);
}

int main()
{
	const bool show_example_eccentric_anomaly_graph = false;
	if(show_example_eccentric_anomaly_graph)
		simple_graph(0, 2*PI);
	
	true_anomaly_over_time();
	
	exit(EXIT_SUCCESS);
}
