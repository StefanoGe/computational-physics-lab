// GenericUtilities.c

#ifndef SG_GENERIC_UTILITIES_C
#define SG_GENERIC_UTILITIES_C

#include <stdlib.h>
#include <stdio.h>
#include "comp_physics.h"

#define UNUSED(expr) do { (void)(expr); } while(0)


#define raiseErr( msg, ...) \
	do { \
		fprintf(stderr, "[%s: line %d in %s]: " msg, __FILE__, __LINE__, __func__ \
		__VA_OPT__(,) __VA_ARGS__); \
		exit(EXIT_FAILURE); \
	} while (0)

/*
#define __OLD_raiseErr(msg, ...) \
    do { \
        fprintf(stderr, msg __VA_OPT__(,) __VA_ARGS__); \
        exit(EXIT_FAILURE); \
    } while(0)
*/

#define eprint(a, ...) fprintf( stderr, a "\n" __VA_OPT__(,) __VA_ARGS__ )

FILE * openFile( char * fileName, char * mode )
{
	FILE * file;
	file = fopen(fileName, mode);
	if (file == NULL) 
		raiseErr( "File %s cannot be opened\n", fileName );
		
    return file;
}

/*
void raiseErr( char * errorMessage )
{
	printf("%s", errorMessage);
	putchar('\n');
	exit(EXIT_FAILURE);
}
*/

#endif
