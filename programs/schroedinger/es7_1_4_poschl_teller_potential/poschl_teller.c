#include "schroedinger.h"
#include "explot.h"

#define MAX_STR_LEN 200

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
	StaticPotential V,void *pot_params, int frameskip, int nframes, double dt,
	char *add_name)
{
	Matrix sol_abs=matc_cabs(sol);
	Matrix sol_real=matc_real_part(sol);
	Matrix sol_imag=matc_imag_part(sol);
	Array xgrid = arr_linspace(x0,xf,nxpoints);
	
	ArrayC potential = arrc_map_par_from_arr(&xgrid,V,pot_params);
	
	Array potential_real = arrc_real(&potential);
	
	AnimSpec cabs_psi = eplot_init_animseries(xgrid.data,nxpoints,sol_abs.rows,
			ntsteps+1,frameskip,"l","|Ψ|^2",NULL,false);
			
	AnimSpec real_psi = eplot_init_animseries(xgrid.data,nxpoints,sol_real.rows,
			ntsteps+1,frameskip,"l","Re(Ψ)",NULL,false);
	
	AnimSpec imag_psi = eplot_init_animseries(xgrid.data,nxpoints,sol_imag.rows,
			ntsteps+1,frameskip,"l","Im(Ψ)",NULL,false);
			
	AnimSpec potential_series = eplot_init_animseries(xgrid.data,nxpoints,
			&(potential_real.data),
			ntsteps+1,frameskip,"l","V",NULL,true);
	
	AnimSpec series[]={cabs_psi,real_psi,imag_psi,potential_series};
	
	char output_name[MAX_STR_LEN];
	snprintf(output_name,MAX_STR_LEN,"poschl_teller_gaussian_%s",add_name);
	
	Global gb_setting = eplot_init_global("Pöschl-Teller potential scattering",
		"x","",NULL,
		output_name,output_name,NULL,
		"set yrange [-0.5:0.5]\n"
			"set xrange[-50:50]");
	
	char set_title[MAX_STR_LEN];
	snprintf(set_title,MAX_STR_LEN,"set title sprintf(\"Gaussian packet "
	"against Pöschl-Teller potential t=%%4.2f\",i*%lf)\n",dt*frameskip);
	
	eplot_anim(series,4,&gb_setting,nframes,
		set_title,
		6);
	
	Matrix *mat_to_free[]={&sol_abs,&sol_imag,&sol_real};
	mat_free_many(mat_to_free,3);
	arr_free(&xgrid);
}

void test_poschl_teller(double omega, double V0, double tf, int ntsteps,
	double xm, double sigma,double k0,bool update_animation,int frameskip,
	int nframes, char *add_name)
{
	double gaussian_params[]={xm,sigma,k0};
	
	double pt_params[]={V0,omega};
	
	const double x0=-100;
	const double xf=100;
	const double t0=0;
	
	const int nxpoints = 2000;
	
	const double dx=(xf-x0)/(nxpoints-1);
	const double dt=(tf-t0)/ntsteps;
	
	MatrixC sol = {0};
	Array time = {0};
	
	schrod_rk4(schrod_gaussian,gaussian_params,schrod_poschl_teller,
		pt_params,&sol,&time,ntsteps,nxpoints,t0,tf,x0,xf);
		
	reflection_and_transmission(&sol,30,dt,-4,4,dx,x0,xf,nxpoints);
		
	if(update_animation)
		animation(&sol,x0,xf,nxpoints,ntsteps,schrod_poschl_teller,
					pt_params,frameskip,nframes,dt,add_name);
}					

int main()
{

//	test_poschl_teller(1,2,40,-20,2,2,false);
//	test_poschl_teller(0.5,-2,80,8000,-20,4,1,true,40,200,"negative");
	const double omega=0.5;
	const double l=4;
	test_poschl_teller(0.5,-l*(l+1)/2/omega/omega,80,8000,-20,4,1,true,40,200,"special");
	
	exit(EXIT_SUCCESS);
}
