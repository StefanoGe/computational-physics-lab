#include "comp_physics.h"

typedef struct {
	char *title;
	char *xlabel;
	char *ylabel;
	char *logscale;
	char *build_name;
	char *output_name;
}Global;

typedef struct{
	char *file;
	char *style;
	char *label;
	char *color;
} DatasetDesc;

bool eplot_set_common(FILE *gp);

Global null_global_settings();

void eplot_set_global(const Global *settings, FILE *gp);

void eplot_2carr(double *xaxis, double *yaxis, int size, const Global *gb_settings,
					const DatasetDesc *data_settings);
