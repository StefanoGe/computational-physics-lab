#include "comp_physics.h"
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

void plot_phase_space(const Matrix *sol, int nsteps)
{
	// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "Lotka-Volterra phase-space");
	sprintf(build_name, "lotkav_phase_sp");
	sprintf(output_name, "lotkav_phase_sp");

	Global gb_settings = {
		.title=title,
		.xlabel="Preys",
		.ylabel="Predators",
		.logscale=nullptr,
		.build_name=build_name,
		.output_name=output_name
	};

	SeriesSpec curve={
		.x=sol->rows[0],
		.y=sol->rows[1],
		.size=nsteps+1,
		.style="l lw 0.5",
		.label="[y,z](0)=[1,0.01]",
		.color=nullptr
	};

	SeriesSpec series[]={curve};
	eplot_multi(series, 1, &gb_settings);
}

void plot_time_evolution(const Array *time, const Matrix *sol)
{
	// Plot setup --- global settings
	char title[1000];
	char build_name[100];
	char output_name[100];
	sprintf(title, "Lotka-Volterra time evolution");
	sprintf(build_name, "lotkav_time_ev");
	sprintf(output_name, "lotkav_time_ev");

	Global gb_settings = {
		.title=title,
		.xlabel="Time",
		.ylabel="Count",
		.logscale=nullptr,
		.build_name=build_name,
		.output_name=output_name
	};

	SeriesSpec preys={
		.x=time->data,
		.y=sol->rows[0],
		.size=time->size,
		.style="l",
		.label="Preys",
		.color=nullptr
	};
	
	SeriesSpec predators={
		.x=time->data,
		.y=sol->rows[1],
		.size=time->size,
		.style="l",
		.label="Predators",
		.color=nullptr
	};

	SeriesSpec series[]={preys, predators};
	eplot_multi(series, 2, &gb_settings);
	
	
}

void study_lv_system(double t0, double tf, double y0, double z0, double alpha0, double beta0, int nsteps)
{
	const double in_value[]={y0,z0};
	Matrix sol={0};
	Array time={0};
	const Derivative eqs[]={eq1, eq2};
	double params[]={alpha0, beta0};
	void *both_params[]={(void*)params, (void*)params};
	
	ode_rk4( t0, tf, in_value, eqs, both_params, 2, nsteps, &time, &sol );
	
	plot_phase_space(&sol, nsteps);
	
	plot_time_evolution(&time, &sol);
	
	Array *arr_to_free[]={&time};
	arr_free_many(arr_to_free, 1);
	Matrix *mat_to_free[]={&sol};
	mat_free_many(mat_to_free, 1);
}

int main()
{
	study_lv_system(0, 60, 1, 0.01, 0.1, 0.25, 1000);
	
	return 0;
}
