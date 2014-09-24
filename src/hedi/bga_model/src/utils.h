void
save_wave_file(char *fname, double *data, int size);

double         *
load_wave_file(char *fname, int *ns);


double
getABSmax(double *data, int size);
void
scaleData(double scale, double *data, int size);
void normalize(double *data,int size);
double
CHI2(double *a, double *b, int size);