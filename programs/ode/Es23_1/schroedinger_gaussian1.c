// Nota: 
// algoritmo adattivo: confronto le soluzioni di una eq diff con la differenza
// fra due metodi, di cui uno meno  euno più preciso. Tengo la soluzione soltanto
// se differenzasarà minore o uguale a una certa tolleranza.
// è un'operazione che faccio a ogni step, e a ogni step riaggiusto lo 
// stepsize
#include "matrixc.h"
#include "arrayc.h"
#include "array.h"
#include <math.h>
#include <string.h>
#include "explot.h"

#define MIN_XPOINTS 100

double nyquist_gaussian(double sigma,double k0,double eps){
	return PI/( k0 + 1/sqrt(2)/sigma*sqrt( log( sigma/eps * sqrt(2/PI) ) ) )/10;}

void init_output(double t0,double tf,int ntimes, int nxsteps,
					Array *t_output, MatrixC *f_output)
{
	if(t_output!=NULL)
	{
		arr_free(t_output);
		*t_output=arr_linspace(t0,tf,ntimes);
	}

	if(f_output->nrows!=ntimes || f_output->ncols!=nxsteps)
		matc_init(f_output, ntimes, nxsteps);
}

// Nota:
// Riassunto:
// Adesso iniziamo solo con V=0
// Procedimento:
// Discretizzo l'hamiltoniana, producendo la matrice
// Applico ripetutamente rk e a ogni step temporale calcolo poi lo step dopo.
// 

void disc_h(MatrixC *M, double xstep_size)
{
	if(!matc_is_squared(M))
		raiseErr("matrix was not initialized correctly");
	
	const int nxpoints = M->nrows;
	
	const double complex a=I/24/xstep_size/xstep_size;
	const double complex g=-I*5/4/xstep_size/xstep_size;
	matc_setv(M, 0);
	for(int i=0; i<nxpoints; i++)
		MATP(M,i,i)=g;
	for(int i=0; i<nxpoints-1; i++)
	{
		MATP(M,i,i+1)=16*a;
		MATP(M,i+1,i)=16*a;
	}
	for(int i=0; i<nxpoints-2; i++)
	{
		MATP(M,i,i+2)=-a;
		MATP(M,i+2,i)=-a;
	}
}

typedef double complex (*WavePacketGen)(double, void*);



double complex gaussian(double x, void *params)
{
	const double xm=((double*)params)[0];
	const double sigma=((double*)params)[1];
	const double k0=((double*)params)[2];
	return exp(-(x-xm)*(x-xm)/4/sigma/sigma)/pow(2*PI*sigma*sigma,0.25)*
		cexp(I*k0*x);
}

void schrod_rk4_core(const MatrixC *M,int nxpoints, int ntsteps,double dt,MatrixC *sol)
{
	ArrayC k1=arrc_new(nxpoints);
	ArrayC k2=arrc_new(nxpoints);
	ArrayC k3=arrc_new(nxpoints);
	ArrayC k4=arrc_new(nxpoints);
	ArrayC temp=arrc_new(nxpoints);
	// apply recursively rk4
	for(int t=0; t<ntsteps; t++)
	{
		// Make k1
		matc_carrc_mult_banded(M, sol->rows[t], &k1,2 );
		// Make k2
		arrc_cp_from_carr(sol->rows[t], nxpoints, &temp);
		arrc_axpy(&temp, dt/2, &k1);
		matc_carrc_mult_banded(M, temp.data, &k2,2);
		// Make k3
		arrc_cp_from_carr(sol->rows[t], nxpoints, &temp);
		arrc_axpy(&temp, dt/2, &k2);
		matc_carrc_mult_banded(M, temp.data, &k3,2);
		// Make k4
		arrc_cp_from_carr(sol->rows[t], nxpoints, &temp);
		arrc_axpy(&temp, dt, &k3);
		matc_carrc_mult_banded(M, temp.data, &k4,2);
		// Update
		for(int j=0; j<nxpoints; j++)
			MATP(sol,t+1,j)=MATP(sol,t,j)+dt/6*(ARR(k1,j)+2*ARR(k2,j)
						+2*ARR(k3,j)+ARR(k4,j));
	}
	
	ArrayC *arrc_to_free[]={&k1,&k2,&k3,&k4,&temp};
	arrc_free_many(arrc_to_free, 5);
	
}

void schrod( WavePacketGen wavepacket, void *wp_params, MatrixC *sol,
			Array *time, int ntsteps, int nxpoints, double t0, double tf,
			double x0, double xf )
{
	const double dx=(xf-x0)/(nxpoints-1);
	const double dt=(tf-t0)/ntsteps;
	const int ntimes=ntsteps+1;
	
	// Allocate output data structs
	init_output(t0,tf,ntimes, nxpoints, time, sol);

	// Discretize hamiltonian. M represents -iH
	MatrixC M=matc_new(nxpoints, nxpoints);
	disc_h(&M, dx);
	
	// setup wavepacket on sol[0]
	for(int i=0; i<nxpoints; i++)
		MATP(sol,0,i)=wavepacket(x0+dx*i, wp_params);
	
	schrod_rk4_core(&M,nxpoints,ntsteps,dt,sol);

	matc_free(&M);
}

void schrod_rk4_gaussian_nyquist(double sigma, double xm, double k0, double eps,
			MatrixC *sol,
			Array *time, int ntsteps, double t0, double tf,
			double x0, double xf)
{
	double wp_params[]={xm,sigma,k0};
	const double xwidth = fabs(xf-x0);
	const int nxpoints = MAX( (int)ceil( nyquist_gaussian(sigma,k0,eps)/xwidth ),
		MIN_XPOINTS );
	schrod(gaussian,(void*)wp_params,sol,time,ntsteps,nxpoints,t0,tf,x0,xf);
}

/*
void test( double xm, double sigma, int nxsteps, int ntsteps,
			double x0, double xf, double t0, double tf )
{
	MatrixC sol;
	Array time;
	double gaussian_params[]={xm, sigma};
	schrod(gaussian, (void*)gaussian_params, &sol, &time, ntsteps, nxsteps,
						t0, tf, x0, xf);
	
	
	
	
	
	matc_free(&sol);
	arr_free(&time);
}
*/
// Remark: we could diagonalize following a routine designed for banded matrices

int main()
{	
	const double xm=0;
	const double sigma=1;
	const double k0=0;
	double gaussian_params[]={xm,sigma,k0};
	const double x0 = -10;
	const double xf =  10;
	const double t0=0;
	const double tf=20;
	
//	const double dx = 0.1;
//	const double dt= 0.01;
	
	const int nxpoints = 200;
	const int ntsteps  = 2000;
	MatrixC sol = {0};
	Array time = {0};
	Array xgrid = arr_linspace(x0,xf,nxpoints);
	
	schrod(gaussian,gaussian_params,&sol,&time,ntsteps,nxpoints,t0,tf,x0,xf);

	Matrix sol_abs=matc_cabs(&sol);
	Matrix sol_real=matc_real_part(&sol);
	Matrix sol_imag=matc_imag_part(&sol);
	
	AnimSpec cabs_psi = eplot_init_animseries(xgrid.data,nxpoints,sol_abs.rows,
			ntsteps+1,10,"l","|Ψ|^2",NULL);
			
	AnimSpec real_psi = eplot_init_animseries(xgrid.data,nxpoints,sol_real.rows,
			ntsteps+1,10,"l","Re(Ψ)",NULL);
	
	AnimSpec imag_psi = eplot_init_animseries(xgrid.data,nxpoints,sol_imag.rows,
			ntsteps+1,10,"l","Im(Ψ)",NULL);
	
	AnimSpec series[]={cabs_psi,real_psi,imag_psi};
	
	Global gb_setting = eplot_init_global("Gaussian wavepacket in a box",
		"x","",NULL,
		"schrodinger_gaussian","schrodinger_gaussian",NULL,"set yrange [-0.7:0.7]");
	
	eplot_anim(series,3,&gb_setting,200,
		"set title sprintf(\"Gaussian wavepacket in a box t=%4.2f\",i*0.01)\n",
		6);
	
	arr_free(&time);
	arr_free(&xgrid);
	matc_free(&sol);
	mat_free(&sol_abs);
	return 0;
}
