#include "ode.h"
#include "explot.h"
#include "array.h"
#include "matrix.h"
#include <math.h>



double eq1([[maybe_unused]] double t, double *u, void *lv_params )
{
	const double y=u[0];
	const double z=u[1];
	const double alpha = ((double*)lv_params)[0];
	const double beta  = ((double*)lv_params)[1];
	return y*(1-alpha*y)-y*z/(1+beta*y);
}

double eq2([[maybe_unused]] double t, double *u, void *lv_params )
{
	const double y=u[0];
	const double z=u[1];

	const double beta  = ((double*)lv_params)[1];
	return -z+y*z/(1+(beta*y));
}

void plot_phase_space(const Matrix *sol, int nsteps, int ncurves,
	const Matrix *closed_orbit, double equilibrium[])
{
	Global gb_settings = eplot_init_global("Predator prey model phase-space",
		"Preys","Predators",NULL,
		"lotkav_phase_sp","lotkav_phase_sp",NULL,
		"set linetype 3 lc rgb \"#B48D30\"");

	SeriesSpec *curves=NULL;
	SAFE_ALLOC(curves,ncurves+2);

	char *labels;
	SAFE_ALLOC(labels,ncurves*1000);
	char (*labels_pointer)[1000]=(char (*)[1000])labels;


	for(int i=0; i<ncurves; i++)
	{
		snprintf(labels_pointer[i], 1000, "[y,z](0)=[%.3g, %.3g]", 
			MATP(sol+i,0,0),MATP(sol+i,1,0));
		curves[i]=eplot_init_series((sol+i)->rows[0],(sol+i)->rows[1],nsteps+1,
			"l lw 0.5",labels_pointer[i],nullptr);
	}
	
	const int closed_points_to_plot = 1000;
	const int first_plot_point_index = closed_orbit->ncols-closed_points_to_plot;
		
	curves[ncurves]=eplot_init_series(closed_orbit->rows[0]+first_plot_point_index,
		closed_orbit->rows[1]+first_plot_point_index, closed_points_to_plot,
		"l lw 2", "Closed orbit", NULL);
		
	curves[ncurves+1]=eplot_init_series(equilibrium,
		equilibrium+1, 1, "p", "Equilibrium point", NULL);

	eplot_multi_options(curves, ncurves+2, &gb_settings, false);
	SAFE_FREE(curves);
}

void plot_time_evolution(const Array *time, const Matrix *sol)
{
	// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "Predator-Prey model time evolution");
	sprintf(build_name, "lotkav_time_ev");
	sprintf(output_name, "lotkav_time_ev");

	Global gb_settings = eplot_init_global(title,"Time","Count",nullptr,build_name,
		output_name,NULL,NULL);

	SeriesSpec preys = eplot_init_series(time->data,sol->rows[0],time->size,
		"l","Preys",nullptr);
	
	SeriesSpec predators=eplot_init_series(time->data,sol->rows[1],time->size,
		"l","Predators",nullptr);

	SeriesSpec series[]={preys, predators};
	eplot_multi_options(series, 2, &gb_settings, false);
}

void study_lv_system(double t0, double tf, const double initial_values[][2], 
	int ncurves, double alpha0, double beta0, int nsteps, 
	double equilibrium[])
{
	Matrix *sol={0};
	SAFE_ALLOC(sol, ncurves);
	Array time={0};
	const Derivative eqs[]={eq1, eq2};
	double params[]={alpha0, beta0};
	void *both_params[]={(void*)params, (void*)params};
	
	for(int i=0; i<ncurves; i++)
		ode_rk4( t0, tf, initial_values[i], eqs, both_params, 2, nsteps, &time, 
			sol+i );
	
	
	Matrix closed_orbit={0};
	Array time_closed={0};
	double final_time_closed=600;
	const double closed_init_values[]={MATP(sol,0,nsteps), MATP(sol,1,nsteps)};
	const int npoints_in_closed_curve = 10000;
	ode_rk4( t0, final_time_closed, closed_init_values, eqs, both_params, 2, 
		npoints_in_closed_curve, &time_closed, &closed_orbit );

	
	plot_phase_space(sol, nsteps, ncurves, &closed_orbit, equilibrium);
	
	plot_time_evolution(&time, sol);
	
	Array *arr_to_free[]={&time};
	arr_free_many(arr_to_free, 1);
	for(int i=0; i<ncurves; i++)
		mat_free(sol+i);
	SAFE_FREE(sol);
}


int main()
{
	double equilibrium[]= {4.0/3.0, 52.0/45.0};
	const double initial_values[][2]={ {1,0.01}, 
		{equilibrium[0]+0.1, equilibrium[1]+0.1}, 
		{5,1} };
	const int npoints = 3;
	study_lv_system(0, 60, initial_values, npoints, 0.1, 0.25, 1000, 
		equilibrium);
	
	return 0;
}
