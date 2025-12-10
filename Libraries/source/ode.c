// ode.c
#include "ode.h"
#include "comp_physics.h"



void ode_feuler( double init_time, double final_time, const Array *init_value, 
				const Derivative *ders, 
				int sys_size, int nsteps, Array *t_output, Matrix *f_output )
{
	if( sys_size!= init_value->size )
		raiseErr("Initial value must be compatible with system's size");

	if(t_output->size!=nsteps+1)
		arr_init(t_output, nsteps+1);
		
	if(f_output->nrows!=nsteps+1 || f_output->ncols!=sys_size)
		mat_init(f_output, nsteps, sys_size);
	
	const double step_size = ( final_time - init_time )/nsteps;
	for(int i=0; i <=nsteps; i++ )
		t_output->data[i]=init_time + step_size*i;
	
	for(int j=0; j<sys_size; j++)
		MATP(f_output,0,j) = init_value->data[j];

	for(int i = 0; i < nsteps; i++ )
	{
		for(int j=0; j<sys_size; j++ )
		{
			MATP(f_output,i+1,j)=MATP(f_output,i,j) 
				+ step_size*ders[j](t_output->data[i], mat_getr_ptr(f_output,i), NULL);
		}
	}
}
