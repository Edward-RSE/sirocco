/* recipes.c */
double num_int (double (*func) (double, void *), double a, double b, double eps);
double zero_find (double (*func) (double, void *), double x1, double x2, double tol, int *ierr);
double find_function_minimum (double a, double m, double b, double (*func) (double, void *), double tol, double *xmin);
int fraction (double value, double array[], int npts, int *ival, double *f, int mode);
int linterp (double x, double xarray[], double yarray[], int xdim, double *y, int mode);
/* random.c */
int randvec (double a[], double r);
int randvcos (double lmn[], double north[]);
double vcos (double x, void *params);
int randvdipole (double lmn[], double north[]);
double vdipole (double cos_theta, void *params);
int init_rand (int seed);
void init_rng_directory (char *root, int rank);
void save_gsl_rng_state (void);
void reload_gsl_rng_state (void);
double random_number (double min, double max);
/* cdf.c */
int cdf_gen_from_func (CdfPtr cdf, double (*func) (double, void *), double xmin, double xmax, int njumps, double jump[]);
double gen_array_from_func (double (*func) (double, void *), double xmin, double xmax, int pdfsteps);
int cdf_gen_from_array (CdfPtr cdf, double x[], double y[], int n_xy, double xmin, double xmax);
double cdf_get_rand (CdfPtr cdf);
int cdf_limit (CdfPtr cdf, double xmin, double xmax);
double cdf_get_rand_limit (CdfPtr cdf);
int cdf_to_file (CdfPtr cdf, char comment[]);
int cdf_inputs_to_file (double x[], double y[], int n_xy, double xmin, double xmax, char filename[]);
int cdf_check (CdfPtr cdf);
int calc_cdf_gradient (CdfPtr cdf);
int cdf_array_fixup (double *x, double *y, int n_xy);
/* vvector.c */
double dot (double a[], double b[]);
double length (double a[]);
int renorm (double a[], double scalar);
int rescale (double a[], double scalar, double b[]);
int cross (double a[], double b[], double c[]);
int vmove (double u[], double lmn[], double s, double result[]);
int vsub (double u[], double v[], double result[]);
int vadd (double u[], double v[], double result[]);
int stuff_v (double vin[], double vout[]);
double dot_tensor_vec (double tensor[3][3], double vin[3], double vout[3]);
int project_from_xyz_cyl (double a[], double b[], double result[]);
int project_from_cyl_xyz (double a[], double b[], double result[]);
int create_basis (double u[], double v[], struct basis *basis_new);
int project_from (struct basis *basis_from, double v_in[], double v_out[]);
int project_to (struct basis *basis_from, double v_in[], double v_out[]);
int reorient (struct basis *basis_from, struct basis *basis_to, double v_from[], double v_to[]);
