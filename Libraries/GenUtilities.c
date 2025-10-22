// GenericUtilities.c

#ifndef SG_GENERIC_UTILITIES_C
#define SG_GENERIC_UTILITIES_C

#include <stdlib.h>
#include <stdio.h>

#define raiseErr(msg, ...) \
    do { \
        fprintf(stderr, msg __VA_OPT__(,) __VA_ARGS__); \
        exit(EXIT_FAILURE); \
    } while(0)

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
