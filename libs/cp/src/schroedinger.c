#include "schroedinger.h"
#include <math.h>
#include "linearsys.h"

#define MIN_XPOINTS 100
//#define DEBUG_RK4 

double complex schrod_gaussian(double x, void *params)
{
	const double xm=((double*)params)[0];
	const double sigma=((double*)params)[1];
	const double k0=((double*)params)[2];
	return exp(-(x-xm)*(x-xm)/4/sigma/sigma)/pow(2*PI*sigma*sigma,0.25)*
		cexp(I*k0*x);
}

double complex schrod_finite_barrier(double x, void *params)
{
	const double V0 = ((double*)params)[0];
	const double x1 = ((double*)params)[1];
	const double x2 = ((double*)params)[2];
	return ( x1<x && x<x2 ? V0 : 0);
}

double complex schrod_poschl_teller(double x, void *params)
{
	const double V0 = ((double*)params)[0];
	const double omega = ((double*)params)[1];
	return V0/pow(cosh(x/omega),2);
}

double nyquist_gaussian(double sigma,double k0,double eps){
	return PI/( k0 + 1/sqrt(2)/sigma*sqrt( log( sigma/eps * sqrt(2/PI) ) ) )/10;}

void schrod_init_output(double t0,double tf,int ntimes, int nxsteps,
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

void schrod_disc_h(MatrixC *M, double x0,double dx, StaticPotential V,
		void *potential_params)
{
	if(!matc_is_squared(M))
		raiseErr("matrix was not initialized correctly");
	
	const int nxpoints = M->nrows;
	
	const double complex a=I/24/dx/dx;
	const double complex g_base=5/4/dx/dx;
	matc_setv(M, 0);
	for(int i=0; i<nxpoints; i++)
		MATP(M,i,i)=-I*(g_base+(V ? V(x0+i*dx,potential_params) : 0));
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

void schrod_disc_free_h(MatrixC *M, double dx)
{
	schrod_disc_h(M,0,dx,NULL,NULL);
}

void schrod_disc_h_improved(MatrixC *P, double x0,double dx, StaticPotential V,
		void *potential_params)
{
	if(!matc_is_squared(P))
		raiseErr("matrix was not initialized correctly");
	
	MatrixC M=matc_new(P->nrows,P->ncols);
	schrod_disc_h(&M,x0,dx,V,potential_params);
	
	
	
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
		#ifdef DEBUG_RK4
			for(int i=0; i<nxpoints; i++)
				printf("%.2lf   %lfe%+lfi\n",i*dt,
					cabs(MATP(sol,t,i)),carg(MATP(sol,t,i)));
			putchar('\n');
			mypause();
		#endif
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

void schrod_rk4_inf_well( WavePacketGen wavepacket, void *wp_params, 
			MatrixC *sol,
			Array *time, int ntsteps, int nxpoints, double t0, double tf,
			double x0, double xf)
{
	const double dx=(xf-x0)/(nxpoints-1);
	const double dt=(tf-t0)/ntsteps;
	const int ntimes=ntsteps+1;
	
	// Allocate output data structs
	schrod_init_output(t0,tf,ntimes, nxpoints, time, sol);

	// Discretize hamiltonian. M represents -iH
	MatrixC M=matc_new(nxpoints, nxpoints);
	schrod_disc_free_h(&M, dx);
	
	// setup wavepacket on sol[0]
	for(int i=0; i<nxpoints; i++)
		MATP(sol,0,i)=wavepacket(x0+dx*i, wp_params);
	
	schrod_rk4_core(&M,nxpoints,ntsteps,dt,sol);

	matc_free(&M);
}

void schrod_rk4( WavePacketGen wavepacket, void *wp_params, 
			StaticPotential V, void *v_params,
			MatrixC *sol,
			Array *time, int ntsteps, int nxpoints, double t0, double tf,
			double x0, double xf )
{
	const double dx=(xf-x0)/(nxpoints-1);
	const double dt=(tf-t0)/ntsteps;
	const int ntimes=ntsteps+1;
	
	// Allocate output data structs
	schrod_init_output(t0,tf,ntimes, nxpoints, time, sol);
	
	// Discretize hamiltonian. M represents -iH
	MatrixC M=matc_new(nxpoints, nxpoints);
	schrod_disc_h(&M, x0, dx, V, v_params);
	
	// setup wavepacket on sol[0]
	for(int i=0; i<nxpoints; i++)
		MATP(sol,0,i)=wavepacket(x0+dx*i, wp_params);
	
	schrod_rk4_core(&M,nxpoints,ntsteps,dt,sol);

	matc_free(&M);
}


void schrod_tise_disc_h(Matrix *H, double x0,double dx, RealPotential V,
		void *potential_params)
{
	if(!mat_is_squared(H))
		raiseErr("matrix was not initialized correctly");
	
	const int nxpoints = H->nrows;
	
	const double b=-1.0/24/dx/dx;
	const double a=5.0/4/dx/dx;
	mat_setv(H, 0);
	for(int i=0; i<nxpoints; i++)
		MATP(H,i,i)=a+(V ? V(x0+i*dx,potential_params) : 0);
	for(int i=0; i<nxpoints-1; i++)
	{
		MATP(H,i,i+1)=16*b;
		MATP(H,i+1,i)=16*b;
	}
	for(int i=0; i<nxpoints-2; i++)
	{
		MATP(H,i,i+2)=-b;
		MATP(H,i+2,i)=-b;
	}
}

void schrod_tise_pure_qr( RealPotential V, void *v_params,
			Array *energies, Matrix *efuncs, int nxpoints,
			double x0, double xf, int n_iter)
{
	const double dx=(xf-x0)/(nxpoints-1);

	// Discretize hamiltonian
	Matrix H=mat_new(nxpoints, nxpoints);
	
	schrod_tise_disc_h(&H, x0, dx, V, v_params);
	
//	mat_print_stdout(&H,"%.2e",true);
	
	linst_pure_qr_iter(&H,energies,efuncs,n_iter, true);

	mat_free(&H);
}

void schrod_tise_shifted_qr (RealPotential V, void *v_params,
			Array *energies, Matrix *efuncs, int nxpoints,
			double x0, double xf, int n_iter)
{
	const double dx=(xf-x0)/(nxpoints-1);

	// Discretize hamiltonian
	Matrix H=mat_new(nxpoints, nxpoints);
	
	schrod_tise_disc_h(&H, x0, dx, V, v_params);
	
//	mat_print_stdout(&H,"%.2e",true);
	
	linst_shifted_qr(&H, energies, efuncs, n_iter, true);

	mat_free(&H);
}

void schrod_tise_shifted_qr_optimized (RealPotential V, void *v_params,
			Array *energies, Matrix *efuncs, int nxpoints,
			double x0, double xf, int n_iter)
{
	const double dx=(xf-x0)/(nxpoints-1);

	// Discretize hamiltonian
	Matrix H=mat_new(nxpoints, nxpoints);
	
	schrod_tise_disc_h(&H, x0, dx, V, v_params);
	
//	mat_print_stdout(&H,"%.2e",true);
	
	linst_shifted_qr_optimized_pentadiagonal(&H, energies, efuncs, n_iter, true);

	mat_free(&H);
}

void schrod_rk4_gaussian_inf_well_nyquist(double sigma, double xm, 
			double k0, double eps,
			MatrixC *sol,
			Array *time, int ntsteps, double t0, double tf,
			double x0, double xf)
{
	double wp_params[]={xm,sigma,k0};
	const double xwidth = fabs(xf-x0);
	const int nxpoints = MAX( (int)ceil( nyquist_gaussian(sigma,k0,eps)/xwidth ),
		MIN_XPOINTS );
	schrod_rk4_inf_well(schrod_gaussian,(void*)wp_params,sol,time,ntsteps,nxpoints,t0,tf,x0,xf);
}

// An implementaion of the trapezoidal rule which is good for our purposes
// since the wavefunction goes on beyond the boundaries anyway

double schrod_wavefunction_norm(complex double *packet,int nxpoints,
				double xstep)
{
	double sum=0;
	for(int i=0; i<nxpoints; i++)
		sum+=cabs(packet[i])*cabs(packet[i]);
	return sum*xstep;
}

ArrayC arrc_map_par_from_arr(const Array *arr, WavePacketGen f,void *params)
{
	ArrayC new = arrc_new(arr->size);
	for(int i=0; i<arr->size; i++)
		ARR(new,i)=f(ARRP(arr,i),params);
	return new;
}


