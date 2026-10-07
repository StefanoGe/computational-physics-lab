#include "schroedinger.h"
#include "explot.h"
#include <math.h>
#include "rootf.h"
#include "linearsys.h"

#define SHOW_ROOT_HELPER false
#define REDO_CALCULATION_PURE false
#define REDO_CALCULATION_SHIFTED true

#define DOMAIN_SIZE 200

typedef enum { 
	PURE,
	SHIFTED
	} QRVariant;

double finite_well(double x, void *params)
{
	double a=((double*)params)[0];
	double V0=((double*)params)[1];
	return (-a<=x && x<=a ? -V0 : 0);
}

double even_E_eq(double E, void *params)
{
	const double a  = ((double*)params)[0];
	const double V0 = ((double*)params)[1];
	const double k = sqrt(2*(E+V0));
	const double pho = sqrt(-2*E);
	return tan(k*a)-pho/k;
}

double odd_E_eq(double E, void *params)
{
	const double a  = ((double*)params)[0];
	const double V0 = ((double*)params)[1];
	const double k = sqrt(2*(E+V0));
	const double pho = sqrt(-2*E);
	return 1.0/tan(k*a)+pho/k;
}

int find_bound_states(Array *energies, int *bound_states_indices)
{
	int nbound=0;
	for(int i=0; i<energies->size; i++)
	{
		if(ARRP(energies,i)<0)
		{
			bound_states_indices[nbound]=i;
			nbound++;
		}
	}
	return nbound;
}

void plot_bound_states(int nbound, int *bound_states_indices, 
		const Matrix *efuncs, double x0, double xf, ParamFuncPtr V,
		void *potential_params, int nxpoints, double dx, const Array *energies,
		const char *mode_name)
{
	SeriesSpec *series;
	SAFE_ALLOC(series,nbound+1);
	
	Array domain=arr_linspace(x0,xf,nxpoints);
	Array potential=arr_map_par(&domain,V,potential_params);
	
	series[0]=eplot_init_series(domain.data,potential.data,potential.size,
			"l","V(x)",NULL);
			
	char (*labels)[50];
	SAFE_ALLOC(labels,nbound);
	
	Array *bound_squared;
	SAFE_ALLOC(bound_squared,nbound);
	
	for(int i=0; i<nbound; i++)
	{
		arr_null(bound_squared+i);
		arr_init(bound_squared+i,nxpoints);
		for(int j=0; j<nxpoints; j++)
		{
			ARRP(bound_squared+i,j)=pow(MATP(efuncs,j,bound_states_indices[i]),2)/dx;
		}
		
		snprintf(labels[i], sizeof labels[i],
			"E = %lf",ARRP(energies,bound_states_indices[i]));
		
		series[i+1]=eplot_init_series(domain.data, bound_squared[i].data,
			nxpoints,"l",labels[i],NULL);
	}

	char title[100];
	snprintf(title, sizeof title, "Bound states of finite well - %s QR algorithm", 
		mode_name );

	Global gb_settings = eplot_init_global(title,"x","|Ψ|^2",NULL,"tise_bound_states",
			"tise_bound_states","","set yrange [-0.5:1]\n"
								   "set xrange [-2:2]");

	eplot_multi(series,nbound+1,&gb_settings);
	
	for (int i = 0; i < nbound; i++)
		arr_free(bound_squared + i);
	SAFE_FREE(bound_squared);
	arr_free(&domain);
	arr_free(&potential);
	SAFE_FREE(series);
	
}

bool should_this_mode_redo_calculation(QRVariant mode)
{
	return (REDO_CALCULATION_PURE && mode == PURE ) || 
		(REDO_CALCULATION_SHIFTED && mode == SHIFTED);
}

void get_energies_and_efuncs(Array *energies, Matrix *efuncs,
		int num_iterations, double x0, double xf, int nxpoints,
		RealPotential potential, void *potential_params, QRVariant mode)
{
	char energies_filename[100];
	char eigenfuncs_filename[100];
	sprintf(energies_filename, "data/energies_%s.dat", mode 
			? "shifted" : "pure");
	sprintf(eigenfuncs_filename, "data/eigenfuncs_%s.dat", mode 
			? "shifted" : "pure");
	
	if( should_this_mode_redo_calculation(mode) )
	{
		FILE *energies_dat=openFile(energies_filename,  "w");
		FILE *efuncs_dat=  openFile(eigenfuncs_filename,"w");
		if(mode == PURE)
			schrod_tise_pure_qr(potential,potential_params,energies,efuncs,
				nxpoints,x0,xf, num_iterations);
		else if(mode == SHIFTED)
		{
			schrod_tise_shifted_qr_optimized(potential, potential_params, energies, efuncs,
				nxpoints, x0, xf, num_iterations);
		}
		
		mat_print(efuncs,efuncs_dat,"% .16lf ","\n");
		arr_print(energies_dat,energies,"% .16lf\n");
		fclose(energies_dat);
		fclose(efuncs_dat);
	}
	else
	{
		*efuncs=mat_new_from_file(eigenfuncs_filename, nxpoints, nxpoints);
		*energies=arr_new_from_file(energies_filename, nxpoints);
	}
}

void exact_energies_graph(double E1, double E2,void *params)
{
	Array domain = arr_linspace(E1,E2,DOMAIN_SIZE);
	Array even_E = arr_map_par(&domain,even_E_eq,params);
	Array odd_E  = arr_map_par(&domain,odd_E_eq, params);

	arr_print_inline(&even_E,"%.2e ",true,true);
	arr_print_inline(&odd_E,"%.2e ",true,true);

	Global gb_settings = eplot_init_global(
		"Helper graph for finding energies roots","E","f(E)",NULL,
		"roots_helper_graph",NULL,NULL,"set yrange [-5:5]");
		
	SeriesSpec even = eplot_init_series(domain.data,even_E.data,domain.size,"l",
		"Even",NULL);
	
	SeriesSpec odd = eplot_init_series(domain.data,odd_E.data,domain.size,"l",
		"Odd",NULL);
		
	SeriesSpec series[] = {even,odd};
	
	eplot_multi(series,2,&gb_settings);
	
	Array *arr_to_free[] = {&domain,&even_E,&odd_E};
	arr_free_many(arr_to_free,3);
}

Array find_exact_energies(void *well_params)
{
	Array exact_energies = arr_new(3);
	Array debug_1 = {0};
	ARR(exact_energies,0) = root_bis(even_E_eq,well_params,-9.9,-8,1e-10,&debug_1);
	ARR(exact_energies,1) = root_bis(odd_E_eq,well_params,-7,-6,1e-10,NULL);
	ARR(exact_energies,2) = root_bis(even_E_eq,well_params,-4,-3,1e-10,NULL);
	
//	arr_print_inline(&exact_energies,"%lf ",true,true);
	
//	arr_print_inline(&debug_1,"%lf ",true,true);

	
	return exact_energies;
}

void compare_numerical_and_exact_energies(const Array *energies, int nbound,
		const int bound_states_indices[], const Array *exact_energies)
{
	puts("\n-------------------------------------------\n"
		 ANSI_GREEN"Comparison between energies found through QR eigenvalue algorithm "
		 "and exact solution\n\n"
		 ANSI_YELLOW"Exact | Numerical | Difference"ANSI_RESET);
	for(int i=0; i<nbound; i++)
	{
		printf("%.2lf | %9.2lf | %.2e\n",ARRP(exact_energies,i),
			ARRP(energies,bound_states_indices[i]),ARRP(exact_energies,i)-
			ARRP(energies,bound_states_indices[i]));
	}
	
	puts("\n-------------------------------------------\n");
}


void test_tise(RealPotential V, void *potential_params, double x0,
	double xf, double dx, int nxpoints, int n_iter, QRVariant mode)
{
	Array energies={0};
	Matrix efuncs={0};
	
	const char *mode_name = ( mode == PURE ? "pure" : "shifted" );
	
	if(should_this_mode_redo_calculation(mode))
		eprint("-----------------------------\n"
			   "Calculation of eigenvalues and eigenvectors with %s qr starts\n"
			   "-----------------------------\n",mode_name);
	
	get_energies_and_efuncs(&energies,&efuncs,n_iter,x0,xf,nxpoints,
		V,(void*)potential_params, mode);
		
	if(should_this_mode_redo_calculation(mode))
		eprint("-----------------------------\n"
			   "Calculation completed\n"
			   "-----------------------------\n");
	
	int	*bound_states_indices;
	SAFE_ALLOC(bound_states_indices,nxpoints);
	
	int nbound = find_bound_states(&energies,bound_states_indices);
	
	plot_bound_states(nbound, bound_states_indices, &efuncs, x0, xf, V,
			potential_params, nxpoints, dx, &energies, mode_name);
	
	Array exact_energies = find_exact_energies(potential_params);
	
	compare_numerical_and_exact_energies(&energies,nbound,bound_states_indices,
		&exact_energies);
	
	arr_free(&exact_energies);

	arr_free(&energies);
	mat_free(&efuncs);
	mypause();	
}

void test_qr_trid(RealPotential V, void *potential_params, double x0,
	double xf, double dx, int nxpoints, int n_iter)
{
	Array energies={0};
	Matrix efuncs={0};
	
	schrod_tise_shifted_qr_optimized(V, potential_params, &energies, &efuncs,
		nxpoints, x0, xf, n_iter);
		
	UNUSED(dx);
	
	arr_free(&energies);
	mat_free(&efuncs);
}

int main()
{	

	const double a=1;
	const double V0=10;
	double well_params[]={a,V0};
	
	if(SHOW_ROOT_HELPER)
		exact_energies_graph(-10,0,(void*)well_params);
	
	const double x0=-5;
	const double xf= 5;
	const double dx=0.01;
	const int nxpoints=(int)round((xf-x0)/dx)+1; // 200
	
//	qr_shifted_test(x0, xf, dx, finite_well, well_params, nxpoints, 1000);
//	test_tise_with_pure_qr(finite_well, well_params, x0, xf, dx, nxpoints, 200);

	test_tise(finite_well, well_params, x0, xf, dx, nxpoints, 1000, SHIFTED);
//	test_tise(finite_well, well_params, x0, xf, dx, nxpoints, 200, PURE);

//	test_qr_trid(finite_well, well_params, x0, xf, dx, 10, 1000);

	

	exit(EXIT_SUCCESS);
}

/*	OUTDATED TEST FOR EIGENVALUES
	Array evec_product=arr_new(nxpoints);
	Array curr_evec=arr_new(nxpoints);
		
	for(int i=0; i<nbound; i++)
	{
		for(int j=0; j<nxpoints; j++)
			ARR(curr_evec,j)=MAT(evecs,j,bound_states_indices[i]);
		mat_vecmult(&H_copy,&curr_evec,&evec_product);
		arr_axpy(&evec_product,-ARR(energies,bound_states_indices[i]),&curr_evec);
		arr_print_inline(&evec_product,"%.2g",true,true);
	}
*/	





