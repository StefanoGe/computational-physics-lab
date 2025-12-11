typedef struct {
	double *data;
	int size;
	bool owns_data;
} Array;

void arr_init( Array *arr, int size );

void arr_free(Array *arr);

Array arr_asarr(double *data, int size);
