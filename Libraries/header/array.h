typedef struct {
	double *data;
	int size;
} Array;

void arr_init( Array *arr, int size );

void arr_free(Array *arr);
