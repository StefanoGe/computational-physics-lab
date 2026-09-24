#include "schroedinger.h"
#include "genutil.h"
#include "explot.h"

#define UPDATE_ANIMATION false

void reflection_and_transmission(const MatrixC *sol, double tmax,double dt,
		double x1,double x2,double dx,double x0,double xf,
		int nxpoints)
{
	const int final_index_reflected = (int)((x1-x0)/dx);
	const int first_index_transmitted = (int)((x2-x0)/dx);
	const int time_index = (int)(tmax/dt);
	const double reflected_norm = schrod_wavefunction_norm(
		sol->rows[time_index],final_index_reflected,dx);
	const double transmitted_norm = schrod_wavefunction_norm(
		sol->rows[time_index]+first_index_transmitted,
			nxpoints-first_index_transmitted,dx);
	const double norm_over_whole_interval = schrod_wavefunction_norm(
		sol->rows[time_index],nxpoints,dx);
	printf( "Norm over the whole interval [%.0lf,%.0lf]: %.3lf\n"
			"R = %lf\n"
			"T = %lf\n"
			"R+T = %lf\n",x0,xf,
			norm_over_whole_interval,
			reflected_norm,
			transmitted_norm,
			reflected_norm+transmitted_norm);
}

void animation(const MatrixC *sol,double x0, double xf,int nxpoints,int ntsteps,
	StaticPotential V,void *pot_params)
{
	Matrix sol_abs=matc_cabs(sol);
	Matrix sol_real=matc_real_part(sol);
	Matrix sol_imag=matc_imag_part(sol);
	Array xgrid = arr_linspace(x0,xf,nxpoints);
	
	ArrayC potential = arrc_map_par_from_arr(&xgrid,V,pot_params);
	
	Array potential_real = arrc_real(&potential);
	
	AnimSpec cabs_psi = eplot_init_animseries(xgrid.data,nxpoints,sol_abs.rows,
			ntsteps+1,10,"l","|Ψ|^2",NULL,false);
			
	AnimSpec real_psi = eplot_init_animseries(xgrid.data,nxpoints,sol_real.rows,
			ntsteps+1,10,"l","Re(Ψ)",NULL,false);
	
	AnimSpec imag_psi = eplot_init_animseries(xgrid.data,nxpoints,sol_imag.rows,
			ntsteps+1,10,"l","Im(Ψ)",NULL,false);
			
	AnimSpec potential_series = eplot_init_animseries(xgrid.data,nxpoints,
			&(potential_real.data),
			ntsteps+1,10,"l","V",NULL,true);
	
	AnimSpec series[]={cabs_psi,real_psi,imag_psi,potential_series};
	
	Global gb_setting = eplot_init_global("Finire barrier scattering",
		"x","",NULL,
		"finite_barrier_gaussian","finite_barrier_gaussian",NULL,"set yrange [-0.5:0.5]\n"
			"set xrange[-50:50]");
	
	eplot_anim(series,4,&gb_setting,400,
		"set title sprintf(\"Gaussian packet against finite barrier t=%4.2f\",i*0.1)\n",
		6);
	
	Matrix *mat_to_free[]={&sol_abs,&sol_imag,&sol_real};
	mat_free_many(mat_to_free,3);
	arr_free(&xgrid);
}

int main()
{
	const double xm=-20;
	const double sigma=2;
	const double k0=2;
	double gaussian_params[]={xm,sigma,k0};
	
	const double V0=2;
	const double x1=-1;
	const double x2=1;
	double barrier_params[]={V0,x1,x2};
	
	const double x0=-100;
	const double xf=100;
	const double t0=0;
	const double tf=40;
	
	const double dx=0.1;
	const double dt=0.01;
	
	const int nxpoints = 2000;
	const int ntsteps  = 4000;
	MatrixC sol = {0};
	Array time = {0};
	
	schrod_rk4(schrod_gaussian,gaussian_params,schrod_finite_barrier,
		barrier_params,&sol,&time,ntsteps,nxpoints,t0,tf,x0,xf);
	
	if(UPDATE_ANIMATION)
		animation(&sol,x0,xf,nxpoints,ntsteps,schrod_finite_barrier,
				barrier_params);
				
	reflection_and_transmission(&sol,30,dt,x1,x2,dx,x0,xf,nxpoints);
	
	exit(EXIT_SUCCESS);
}
