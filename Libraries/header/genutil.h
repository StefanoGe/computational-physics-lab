#ifndef GENUTIL_H
#define GENUTIL_H

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#ifndef PI
#define PI 3.1415926535897932384626433
#endif

#define UNUSED(expr) do { (void)(expr); } while(0)

#define raiseErr( msg, ...) \
	do { \
		fprintf(stderr, "[%s: line %d in %s]: " msg, __FILE__, __LINE__, __func__ \
		__VA_OPT__(,) __VA_ARGS__); \
		exit(EXIT_FAILURE); \
	} while (0)

#define eprint(a, ...) fprintf( stderr, a "\n" __VA_OPT__(,) __VA_ARGS__ )

FILE * openFile( const char * fileName, const char * mode );

void mypause ( void );

int factorial(int n);

double pown( double x, int n );

void carr_print_inline(const double *arr, int size, const char *format, 
	bool parentheses, bool linebreak);

#endif
