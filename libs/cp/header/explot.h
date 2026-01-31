#ifndef EXPLOT_H
#define EXPLOT_H

#define GP_INFINITY "{/Symbol oo}"

typedef struct {
	char *title;
	char *xlabel;
	char *ylabel;
	char *logscale;
	char *build_name;
	char *output_name;
	char *key;
	char *additional_set;
}Global;

typedef struct{
	char *style;
	char *label;
	char *color;
} DatasetDesc;

typedef struct{
	double *x;
	double *y;
	int size;
	char *style;
	char *label;
	char *color;
} SeriesSpec;

static const Global NULL_GLOBAL_SETTINGS={nullptr, nullptr, nullptr, nullptr, 
				nullptr, nullptr, nullptr, nullptr};

static const DatasetDesc NULL_DATADESC={nullptr, nullptr, nullptr};

void eplot_set_common(FILE *gp);

Global null_global_settings();

DatasetDesc null_datadesc();

void eplot_set_global( FILE *gp, const Global *settings);

void eplot_2carr(double *xaxis, double *yaxis, int size, const Global *gb_settings,
					const DatasetDesc *data_settings);

void eplot_multi(const SeriesSpec *series, int count, const Global *gb_settings);

#endif
