#include "comp_physics.h"
#include <string.h>

int main()
{
	VectorD domain = vec_range(0, 100, 100);
	VectorD y_values = vec_range(0, 1000, 100);
	tmplot_2vecs( "test", &domain, &y_values );
	
	free_vecD( &domain, &y_values );
	return 0;
}
