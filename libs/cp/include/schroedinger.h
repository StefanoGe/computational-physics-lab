#ifndef SCHROEDINGER_H
#define SCHROEDINGER_H

#include "matrixc.h"
#include "genutil.h"

typedef double complex (*WavePacketGen)(double, void*);

typedef double complex (*StaticPotential)(double, void*);

typedef ParamFuncPtr RealPotential;

void schrod_rk4_inf_well( WavePacketGen wavepacket, void *wp_params, 
			MatrixC *sol,
			Array *time, int ntsteps, int nxpoints, double t0, double tf,
			double x0, double xf );
			
void schrod_rk4( WavePacketGen wavepacket, void *wp_params, 
			StaticPotential V, void *v_params,
			MatrixC *sol,
			Array *time, int ntsteps, int nxpoints, double t0, double tf,
			double x0, double xf );

double complex schrod_gaussian(double x, void *params);

double complex schrod_finite_barrier(double x, void *params);

double complex schrod_poschl_teller(double x, void *params);

double schrod_wavefunction_norm(complex double *packet,int nxpoints,
				double xstep);
				
ArrayC arrc_map_par_from_arr(const Array *arr, WavePacketGen f,void *params);

void schrod_tise_pure_qr( RealPotential V, void *v_params,
			Array *energies, Matrix *efuncs, int nxpoints,
			double x0, double xf, int n_iter);

void schrod_tise_disc_h(Matrix *H, double x0,double dx, RealPotential V,
		void *potential_params);
		

void schrod_tise_shifted_qr( RealPotential V, void *v_params,
			Array *energies, Matrix *efuncs, int nxpoints,
			double x0, double xf, int n_iter);

void schrod_tise_shifted_qr_optimized (RealPotential V, void *v_params,
			Array *energies, Matrix *efuncs, int nxpoints,
			double x0, double xf, int n_iter);

#endif
