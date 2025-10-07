// GenericUtilities.c

FILE * openFile( char * fileName, char * mode )
{
	FILE * file;
	file = fopen(fileName, mode);
	if (file == NULL) {
		printf("The file %s is not opened.", fileName);
    }
    return file;
}
