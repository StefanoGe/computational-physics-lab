#include "comp_physics.h"
#include <string.h>

int main()
{
	char buffer[100];
	config_parser( "test.cfg", buffer );
	printf("%s\n", buffer);
	
	return 0;
}



