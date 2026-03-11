#include <stdlib.h>
#include <stdio.h>
#include "linearsys.h"
#include "genutil.h"

void test_min_x(  )
{
	Matrix test_mat = mat_new_from_file( "test_mat_ls.txt", 3, 2 );
	double test_b []= { 1.0, -5.0, 6.0 };
	Array test_b_arr = arr_asarr(test_b, 3);
	
	Array min_x={0};

	linst_lsqr_lup( &test_mat, &test_b_arr, &min_x );

	printf("x which minimizes required least square system:\n");
	arr_print_inline( &min_x, "%g", true, true );
	
	mat_free(&test_mat);
	arr_free(&min_x);
}

int main()
{
	test_min_x();
	
	return 0;
}





