// GenericUtilities.c

#ifndef SG_GENERIC_UTILITIES_C
#define SG_GENERIC_UTILITIES_C

FILE * openFile( char * fileName, char * mode )
{
	FILE * file;
	file = fopen(fileName, mode);
	if (file == NULL) {
		printf("The file %s is not opened.", fileName);
    }
    return file;
}

void raiseErr( char * errorMessage )
{
	printf("%s", errorMessage);
	putchar('\n');
	exit(EXIT_FAILURE);
}

#endif
