// ode.c
#include "ode.h"
#include "comp_physics.h"

// Common helper functions

void init_output(int sys_size, int nsteps, Array *t_output, Matrix *f_output)
{
	if(t_output->size!=nsteps+1)
		arr_init(t_output, nsteps+1);

	if(f_output->nrows!=nsteps+1 || f_output->ncols!=sys_size)
		mat_init(f_output, nsteps+1, sys_size);
}

void init_starting_values(double step_size, int sys_size, int nsteps, 
					Array *t_output, Matrix *f_output,double init_time,
							const double init_value[])
{
	for(int i=0; i <=nsteps; i++ )
		ARRP(t_output,i)=init_time + step_size*i;

	for(int j=0; j<sys_size; j++)
		MATP(f_output,0,j) = init_value[j];	
}

void **setup_parameters( int sys_size, void **f_parameters )
{
	void **parameters=malloc(sizeof(void*)*sys_size);
	for(int i=0; i<sys_size;i++)
		parameters[i]=(f_parameters)?f_parameters[i]:nullptr;
	return parameters;
}

/*
 * Forward Euler method for solving systems of differential equations that
 * can be written in the form ui'(t)=f(u[],t)
 * 
 * Paramters:
 * init_time      
 * final_time      
 * init_value[]   - value of the functions at time t
 * ders[]         - array of parametrized functions. ders[i](t, u, f_parameters[i]) 
 *                  should return the value of the derivative of the i-th 
 *                  function at time t
 * **f_paramters  - optional parameters of the functions
 * sys_size       - number of equations to be considered, and of functions to be
 *                  computed
 * nsteps         - number of timesteps to divide the time interval equally
 * *t_output      - array to store the initial time + nsteps timesteps 
 * *f_output      - matrix to store the values of the numerical solution.
 *                  f_output.data[i][j] will be the value of the i-th function
 *                  at time indicised by j 
 * 
 */ 

void ode_feuler( double init_time, double final_time, const double init_value[], 
				const Derivative ders[], void **f_parameters,
				int sys_size, int nsteps, Array *t_output, Matrix *f_output )
{
	init_output(sys_size, nsteps, t_output, f_output);
	const double step_size = ( final_time - init_time )/nsteps;
	init_starting_values(step_size,sys_size, nsteps, t_output, f_output,
								init_time, init_value);
	void **parameters=setup_parameters(sys_size, f_parameters);

	for(int t = 0; t < nsteps; t++ )
	{
		for(int j=0; j<sys_size; j++ )
		{
			MATP(f_output,t+1,j)=MATP(f_output,t,j) 
				+ step_size*ders[j](ARRP(t_output,t), f_output->rows[t], parameters);
			if(t<100)
				eprint("component %d at time %d: %lf", j, t+1,MATP(f_output,t+1,j));
		}
	}

	mat_transpose(f_output,f_output);
}




void ode_ie2( double init_time, double final_time, const double init_value[], 
				const Derivative ders[], void **f_parameters,
				int sys_size, int nsteps, Array *t_output, Matrix *f_output )
{
	init_output(sys_size, nsteps, t_output, f_output);
	const double step_size = ( final_time - init_time )/nsteps;
	init_starting_values(step_size,sys_size, nsteps, t_output, f_output,
								init_time, init_value);
	void **parameters=setup_parameters(sys_size, f_parameters);

	Array step1=arr_new(sys_size);
	for(int t = 0; t < nsteps; t++ )
	{
		for(int j=0; j<sys_size; j++)
		{
			ARR(step1,j)=MATP(f_output,t,j) 
				+ step_size/2*ders[j](ARRP(t_output,t), f_output->rows[t], parameters[j]);
		}
		
		const double t_step2=ARRP(t_output,t)+step_size/2;
		
		for(int j=0; j<sys_size; j++ )
		{		
			MATP(f_output,t+1,j)=MATP(f_output,t,j) 
				+ step_size*ders[j](t_step2, step1.data, parameters[j]);
			/*if(t<100)
				eprint("component %d at time %d: %lf", j, t+1,MATP(f_output,t+1,j));*/
		}
	}
	arr_free(&step1);
	free(parameters);
	mat_transpose(f_output,f_output);
}



void ode_rk4( double init_time, double final_time, const double init_value[], 
				const Derivative ders[], void **f_parameters,
				int sys_size, int nsteps, Array *t_output, Matrix *f_output )
{
	init_output(sys_size, nsteps, t_output, f_output);
	const double step_size = ( final_time - init_time )/nsteps;
	init_starting_values(step_size,sys_size, nsteps, t_output, f_output,
								init_time, init_value);
	void **parameters=setup_parameters(sys_size, f_parameters);
	
	Array step1=arr_new(sys_size);
	Array step2=arr_new(sys_size);
	Array step3=arr_new(sys_size);
	
	for(int t = 0; t < nsteps; t++ )
	{
		const double t_step23=ARRP(t_output,t)+step_size/2;
		const double t_step4=ARRP(t_output,t)+step_size;
		// step 1
		for(int j=0; j<sys_size; j++)
			ARR(step1,j)=MATP(f_output,t,j)+step_size/2*ders[j](ARRP(t_output,t), 
								f_output->rows[t], parameters[j]);
		// step 2
		for(int j=0;j<sys_size;j++)
			ARR(step2,j)=MATP(f_output,t,j)+step_size/2*ders[j](t_step23, 
								step1.data, parameters[j]);
		// step 3
		for(int j=0;j<sys_size;j++)
			ARR(step3,j)=MATP(f_output,t,j)+step_size*ders[j](t_step23, 
								step2.data, parameters[j]);
		
		for(int j=0; j<sys_size; j++ )
		{		
			MATP(f_output,t+1,j)=-MATP(f_output,t,j)/3+ARR(step1,j)/3
				+ARR(step2,j)*2/3+ARR(step3,j)/3+
				+step_size*ders[j](t_step4, step3.data, parameters[j])/6;
			/*if(t<100)
				eprint("component %d at time %d: %lf", j, t+1,MATP(f_output,t+1,j));*/
		}
	}
	Array *arr_to_free[]={&step1, &step2, &step3};
	arr_free_many(arr_to_free,3);
	free(parameters);
	mat_transpose(f_output,f_output);
}



