#include "comp_physics.h"
#include <stdlib.h>
#include <stdio.h>



int main()
{
	VectorD points = build_vecD( 4, 0.0, 1.0, 2.0, 3.0 );
	
	std_print_vecD( &points );
	
	free_vecD(&points );
	
	
	exit(EXIT_SUCCESS);
}
