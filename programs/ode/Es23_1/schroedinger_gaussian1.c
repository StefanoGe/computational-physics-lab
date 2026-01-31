// Nota: 
// algoritmo adattivo: confronto le soluzioni di una eq diff con la differenza
// fra due metodi, di cui uno meno  euno più preciso. Tengo la soluzione soltanto
// se differenzasarà minore o uguale a una certa tolleranza.
// è un'operazione che faccio a ogni step, e a ogni step riaggiusto lo 
// stepsize
#include "comp_physics.h"
#include "matrixc.h"
#include "arrayc.h"
#include "array.h"
#include <math.h>
#include <string.h>

void init_output(int ntsteps, int nxsteps,
					Array *t_output, MatrixC *f_output)
{
	if(t_output->size!=ntsteps+1)
		arr_init(t_output, ntsteps+1);

	if(f_output->nrows!=ntsteps+1 || f_output->ncols!=nxsteps)
		matc_init(f_output, ntsteps+1, nxsteps);
}

// Nota:
// Riassunto:
// Adesso iniziamo solo con V=0
// Procedimento:
// Discretizzo l'hamiltoniana, producendo la matrice
// Applico ripetutamente rk e a ogni step temporale calcolo poi lo step dopo.
// 

void disc_h(MatrixC *M, double xstep_size, int nxsteps)
{
	const double complex a=I/24/xstep_size/xstep_size;
	const double complex g=-I*5/4/xstep_size/xstep_size;
	matc_set_value(M, 0);
	for(int i=0; i<nxsteps; i++)
		MATP(M,i,i)=g;
	for(int i=0; i<nxsteps-1; i++){
		MATP(M,i,i+1)=16*a;
		MATP(M,i+1,i)=16*a;
	}
	for(int i=0; i<nxsteps-2; i++){
		MATP(M,i,i+2)=-a;
		MATP(M,i+2,i)=-a;
	}
}

typedef double (*WavePacket)(double, void*);

double gaussian(double x, void *params)
{
	const double xm=((double*)params)[0];
	const double sigma=((double*)params)[1];	
	return exp(-(x-xm)*(x-xm)/4/sigma/sigma)/pow(2*PI*sigma*sigma,0.25);
}


void schrod( WavePacket wavepacket, void *wp_params, MatrixC *sol,
			Array *time, int ntsteps, int nxsteps, double t0, double tf,
			double x0, double xf)
{
	const double dx=(xf-x0)/(nxsteps-1);
	const double dt=(tf-t0)/ntsteps;
	const int ntimes=ntsteps+1;
	
	// Allocate output data structs
	init_output(ntsteps, nxsteps, time, sol);

	// Discretize hamiltonian. M represents -iH
	MatrixC M=matc_new(nxsteps, nxsteps);
	disc_h(&M, dx, nxsteps);
	
	// Fill time
	for(int t=0; t<ntimes; t++)
		ARRP(time,t)=t0+t*dt;
	
	// setup wavepacket in some way, on sol[0]
	for(int i=0; i<nxsteps; i++)
		MATP(sol,0,i)=wavepacket(x0+dx*i, wp_params);
	
	ArrayC k1=arrc_new(nxsteps);
	ArrayC k2=arrc_new(nxsteps);
	ArrayC k3=arrc_new(nxsteps);
	ArrayC k4=arrc_new(nxsteps);
	ArrayC temp=arrc_new(nxsteps);
	// apply recursively rk4
	for(int t=0; t<ntsteps; t++)
	{
		// Make k1
		matc_arrc_mult(M.rows, sol->rows[t], nxsteps, nxsteps, &k1 );
		// Make k2
		arrc_cp(sol->rows[t], nxsteps, &temp);
		arrc_add_scaled(&temp, dt/2, &k1);
		matc_arrc_mult(M.rows, temp.data, nxsteps, nxsteps, &k2);
		// Make k3
		arrc_cp(sol->rows[t], nxsteps, &temp);
		arrc_add_scaled(&temp, dt/2, &k2);
		matc_arrc_mult(M.rows, temp.data, nxsteps, nxsteps, &k3);
		// Make k4
		arrc_cp(sol->rows[t], nxsteps, &temp);
		arrc_add_scaled(&temp, dt, &k3);
		matc_arrc_mult(M.rows, temp.data, nxsteps, nxsteps, &k4);
		// Update
		for(int j=0; j<nxsteps; j++)
			MATP(sol,t+1,j)=MATP(sol,t,j)+dt/6*(ARR(k1,j)+2*ARR(k2,j)
						+2*ARR(k3,j)+ARR(k4,j));
	}
	
	ArrayC *arrc_to_free[]={&k1,&k2,&k3,&k4,&temp};
	arrc_free_many(arrc_to_free, 5);
	matc_free(&M);
}


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

// Remark: we could diagonalize following a routine designed for banded matrices

int main()
{

	
	
	return 0;
}


