#include "matrixc.h"
#include "arrayc.h"
#include "array.h"
#include <math.h>
#include <string.h>
#include "explot.h"
#include "schroedinger.h"

#define WIDER_GRID_WIDTH 1'000

// Question: is the error used in the integration rule used relevant?

void norm_analysis(const MatrixC *sol, int ntsteps, int nxpoints, double dx,
		double dt,double *gaussian_params)
{
	Array wider_grid=arr_linspace(-WIDER_GRID_WIDTH,WIDER_GRID_WIDTH,
		(int)(2*WIDER_GRID_WIDTH/dx));
	
	ArrayC wider_starting_wp = arrc_map_par_from_arr(&wider_grid,schrod_gaussian,
		(void*)gaussian_params);
	
	double full_norm = schrod_wavefunction_norm(wider_starting_wp.data,
		wider_starting_wp.size,dx);
		
	printf( "Without setting the function at 0 outside the well,the norm-1 at "
			"time 0 is:\n %e\n",full_norm-1);
	
	double starting_norm=schrod_wavefunction_norm(sol->rows[0],nxpoints,dx);

	for(int i=0; i<=ntsteps; i++)
		printf("At time t=%.2lf diff with starting norm is %+.4e\n",i*dt,
				schrod_wavefunction_norm(sol->rows[i],nxpoints,dx)-starting_norm);
				
	arr_free(&wider_grid);
	arrc_free(&wider_starting_wp);
}

void animation(const MatrixC *sol,double x0, double xf,int nxpoints,int ntsteps)
{
	Matrix sol_abs=matc_cabs(sol);
	Matrix sol_real=matc_real_part(sol);
	Matrix sol_imag=matc_imag_part(sol);
	Array xgrid = arr_linspace(x0,xf,nxpoints);
	
	AnimSpec cabs_psi = eplot_init_animseries(xgrid.data,nxpoints,sol_abs.rows,
			ntsteps+1,10,"l","|Ψ|^2",NULL,false);
			
	AnimSpec real_psi = eplot_init_animseries(xgrid.data,nxpoints,sol_real.rows,
			ntsteps+1,10,"l","Re(Ψ)",NULL,false);
	
	AnimSpec imag_psi = eplot_init_animseries(xgrid.data,nxpoints,sol_imag.rows,
			ntsteps+1,10,"l","Im(Ψ)",NULL,false);
	
	AnimSpec series[]={cabs_psi,real_psi,imag_psi};
	
	Global gb_setting = eplot_init_global("\"Moving\" gaussian wavepacket in a box",
		"x","",NULL,
		"schrodinger_gaussian2","schrodinger_gaussian2",NULL,"set yrange [-0.5:0.5]");
	
	eplot_anim(series,3,&gb_setting,200,
		"set title sprintf(\"Gaussian wavepacket in a box t=%4.2f\",i*0.1)\n",
		6);
	
	Matrix *mat_to_free[]={&sol_abs,&sol_imag,&sol_real};
	mat_free_many(mat_to_free,3);
	arr_free(&xgrid);
}

int main()
{	
	const double xm=20;
	const double sigma=2;
	const double k0=2;
	double gaussian_params[]={xm,sigma,k0};
	const double x0 = 0;
	const double xf = 100;
	const double t0=0;
	const double tf=20;
	
	const double dx = 0.1;
	const double dt= 0.01;
	
	const int nxpoints = 1000;
	const int ntsteps  = 2000;
	MatrixC sol = {0};
	Array time = {0};
	
	schrod_rk4_inf_well(schrod_gaussian,
		gaussian_params,&sol,&time,ntsteps,nxpoints,t0,tf,x0,xf);

	norm_analysis(&sol,ntsteps,nxpoints,dx,dt,gaussian_params);

	animation(&sol,x0,xf,nxpoints,ntsteps);
	
	arr_free(&time);
	matc_free(&sol);
	return 0;
}

