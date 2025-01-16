#include "hawaii.h"

double pes(double *q) { UNUSED(q); return 0; }
void dpes(double *q, double *dVdq) { UNUSED(q); UNUSED(dVdq); }
    
// WingParams INIT_WP;
    

/*
 * MODEL: Lorentzian function shifted upwards by constant: 
 *             y = C + A /(1 + B^2 x^2)
 */
typedef struct {
    double A;
    double B;
    double C;
} WingParams;
    
typedef struct {
    size_t n;
    double* t;
    double* y;
} WingData;
    
WingParams INIT_WP = {
  .A = 1.0, 
  .B = 2.0, 
  .C = 3.0,
};
    
int wingmodel_f(const gsl_vector* x, void* data, gsl_vector* f)
/*
 * [input]  x    : parameters
 * [input]  data : CF
 * [output] f    : Yi - yi = difference between model and CF values
 *
 * MODEL: Lorentzian function shifted upwards by constant: f = c + a /(1 + b^2 x^2)
 * https://mathworld.wolfram.com/LorentzianFunction.html 
 */ 
{
    WingData *wd = (WingData*) data; 
    size_t n  = wd->n;
    double* t = wd->t;
    double* y = wd->y;

    double A = gsl_vector_get(x, 0);
    double B = gsl_vector_get(x, 1);
    double C = gsl_vector_get(x, 2);

    for (size_t k = 0; k < n; ++k) {
        double Yk = C + A / (1.0 + B * B * t[k] * t[k]);
        gsl_vector_set(f, k, Yk - y[k]);
    }

    return GSL_SUCCESS;
}

int wingmodel_df(const gsl_vector* x, void* data, gsl_matrix * J)
/*
 * [input]  x    : parameters
 * [input]  data : CF
 * [output] J    : Jacobian matrix J(i, j) = dfi / dxj
 *                 where fi = (Yi - yi) -- difference between model and CF values
 *                       and the xj are the parameters  
 *
 * MODEL: Lorentzian function shifted upwards by constant: f = c + a /(1 + b^2 x^2)
 * https://mathworld.wolfram.com/LorentzianFunction.html 
 */ 
{
    size_t n = ((WingData*) data)->n;
    double* t = ((WingData*) data)->t;

    double A = gsl_vector_get(x, 0);
    double B = gsl_vector_get(x, 1);

    for (size_t k = 0; k < n; ++k) {
        gsl_matrix_set(J, k, 0,  1.0 / (1.0 + B * B * t[k] * t[k]));
        gsl_matrix_set(J, k, 1, -2.0 * A * B * t[k] * t[k] / (B * B * t[k] * t[k] + 1) / (B * B * t[k] * t[k] + 1));
        gsl_matrix_set(J, k, 2,  1.0);
    }

    return GSL_SUCCESS;
}

void gsl_multifit_callback(const size_t iter, void* params, const gsl_multifit_nlinear_workspace* w)
{
    UNUSED(params);

    gsl_vector* f = gsl_multifit_nlinear_residual(w);
    gsl_vector* x = gsl_multifit_nlinear_position(w);
    double rcond;

    // compute reciprocal condition number of J(x) 
    gsl_multifit_nlinear_rcond(&rcond, w);

    fprintf(stdout, "    [callback] LM iter %2zu: A = %.12f, B = %.12f, C = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n",
            iter, gsl_vector_get(x, 0), gsl_vector_get(x, 1), gsl_vector_get(x, 2), 1.0 / rcond, gsl_blas_dnrm2(f));
}


void gsl_nonlinear_opt(size_t n, double* x, double* y, WingParams *wing_params)
/*
 * fit wing model to [x, y] data using the Levenberg-Marquardt method
 *
 * based on example from
 * https://www.gnu.org/software/gsl/doc/html/nls.html 
 *
 * WingParams = {A, B, C} : initial parameters provided to LM method
 *                          optimized parameters are written back into array 
 *
 */
{
#define nparams 3 

    printf("n = %zu\n", n);

    WingData d = {
        .n = n, 
        .t = x, 
        .y = y,
    };

    const gsl_multifit_nlinear_type* T = gsl_multifit_nlinear_trust;

    gsl_multifit_nlinear_fdf fdf;
    fdf.f      = wingmodel_f;
    fdf.df     = wingmodel_df; // set to NULL for finite-difference Jacobian
    fdf.fvv    = NULL; // not using geodesic acceleration
    fdf.n      = n;
    fdf.p      = nparams;
    fdf.params = (void*) &d;
    
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
    gsl_multifit_nlinear_workspace* w = gsl_multifit_nlinear_alloc(T, &fdf_params, n, nparams);
    
    double *weights = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; ++i) { 
        weights[i] = 1.0;
    }

    double _wp[nparams] = {wing_params->A, wing_params->B, wing_params->C};
    gsl_vector_view wp = gsl_vector_view_array(_wp, nparams);

    gsl_vector_view wts = gsl_vector_view_array(weights, n);
    gsl_multifit_nlinear_winit(&wp.vector, &wts.vector, &fdf, w); // initialize solver with starting point and weights

    // compute initial cost function
    gsl_vector* f = gsl_multifit_nlinear_residual(w);

    double chisq0;
    gsl_blas_ddot(f, f, &chisq0);

    int status, info;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 1e-8;

    const size_t niter_max = 100;
    status = gsl_multifit_nlinear_driver(niter_max, xtol, gtol, ftol, gsl_multifit_callback, NULL, &info, w);

    // covariance of best fit parameters
    gsl_matrix* J = gsl_multifit_nlinear_jac(w);
    gsl_matrix* covar = gsl_matrix_alloc(3, 3);
    gsl_multifit_nlinear_covar(J, 0.0, covar);

    double chisq;
    gsl_blas_ddot(f, f, &chisq); // final cost

    fprintf(stderr, "    summary from method '%s/%s'\n", gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "    number of iterations: %zu\n", gsl_multifit_nlinear_niter(w));
    fprintf(stderr, "    function evaluations: %zu\n", fdf.nevalf);
    fprintf(stderr, "    Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "    reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "    initial |f(x)| = %f\n", sqrt(chisq0));
    fprintf(stderr, "    final   |f(x)| = %f\n", sqrt(chisq));

    double dof = n - nparams;
    double c = fmax(1.0, sqrt(chisq / dof));
    
    wing_params->A = gsl_vector_get(w->x, 0);
    wing_params->B = gsl_vector_get(w->x, 1);
    wing_params->C = gsl_vector_get(w->x, 2);

    fprintf(stderr, "chisq / dof = %g\n", chisq / dof);
    fprintf(stderr, "A = %.10f +/- %.10f\n", wing_params->A, c * sqrt(gsl_matrix_get(covar, 0, 0)));
    fprintf(stderr, "B = %.10f +/- %.10f\n", wing_params->B, c * sqrt(gsl_matrix_get(covar, 1, 1)));
    fprintf(stderr, "C = %.10f +/- %.10f\n", wing_params->C, c * sqrt(gsl_matrix_get(covar, 2, 2)));

    fprintf(stderr, "status = %s\n", gsl_strerror(status));

    gsl_multifit_nlinear_free(w);
    gsl_matrix_free(covar);

#undef nparams 
}



WingParams fit_baseline(CFnc *cf, size_t EXT_RANGE_MIN)
// EXT_RANGE_MIN :: leftmost value of `tail` range to fit the parameters of the baseline function

// See paper "Simulation of collision-induced absorption spectra based on classical trajectories and
// ab initio potential and induced dipole surfaces. II. CO2-Ar rototranslational band including true 
// dimer contribution"
//
// we rely on the specific structure of the correlation function (see Fig. 7)
// specifically, we expect it to have a secondary maximum from which it decays to some POSITIVE constant
{
    double CFmin  = FLT_MAX; 
    double CFmax  = FLT_MIN; 

    // auxiliary variable to detect that we traversed the first minimum of the CF
    bool wasnegative = false;

    // secondary maximum: its time coordinate [max2time] and CF value [CFmax2]
    double max2time = 0.0; 
    double CFmax2 = FLT_MIN; 

    for (size_t k = 0; k < cf->len; ++k) {
        if (cf->data[k] > CFmax) CFmax = cf->data[k];
        if (cf->data[k] < CFmin) CFmin = cf->data[k];

        if (cf->data[k] < 0.0) wasnegative = true;
        if (wasnegative && (cf->data[k] > CFmax2)) { CFmax2 = cf->data[k]; max2time = cf->t[k]; }
    } 

    printf("INFO: CFmin    = %.10e\n", CFmin); 
    printf("INFO: CFmax    = %.10e\n", CFmax);
    printf("INFO: max2time = %.10e\n", max2time / ATU); 
    printf("INFO: CFmax2   = %.10e\n", CFmax2);
    
    assert(cf->len > EXT_RANGE_MIN);
    
    size_t itotal = cf->len - EXT_RANGE_MIN + 1;
    double *xdat2 = (double*) malloc(itotal * sizeof(double));
    double *ydat2 = (double*) malloc(itotal * sizeof(double));
   
    // selecting values of `wing tail` and also scaling them  
    for (size_t k = EXT_RANGE_MIN, ind = 0; k < cf->len; ++k, ++ind) {
        xdat2[ind] = cf->t[k] / max2time;
        ydat2[ind] = cf->data[k] / CFmax2;
    }

    WingParams wp = {
        .A = INIT_WP.A, 
        .B = INIT_WP.B, 
        .C = fmax(ydat2[itotal-1], INIT_WP.C)
    }; 
    
    printf("Initial parameter values: %.5e %.5e %.5e\n", wp.A, wp.B, wp.C);

    gsl_nonlinear_opt(itotal, xdat2, ydat2, &wp);

    //assert(wp.A > 0);
    //assert(wp.B > 0);
    //assert(wp.C > 0);
    if (wp.C < 0) {
        printf("\n");
        printf("!!! WingParams.C is negative! Deal with it. Continuing...\n\n");
    }

    max2time = max2time / ATU;

    // descale parameters after optimization
    wp.A *= CFmax2;
    wp.B /= max2time;
    wp.C *= CFmax2;

    return wp; 
}

int main()
{
    String_Builder sb = {0};
    CFnc cf = {0};

    const char *filename = "examples/CF-F-300.0.txt"; 
    if (!read_correlation_function(filename, &sb, &cf)) {
        printf("ERROR: could not read the file!\n");
    } 
   
    printf("INFO: Loaded correlation function from %s\n", filename);
    printf("HEADER:\n %s\n\n", sb.items);
    
    for (size_t i = 0; i < cf.len; ++i) {
        cf.t[i] = cf.t[i] * ATU; 
    }

    printf("Number of samples: %zu\n", cf.len);
    printf("t range: %.3lf...%.3lf\n", cf.t[0], cf.t[cf.len - 1]);
    
    size_t EXT_RANGE_MIN = 8192; // far-infrared: 70K-300K
    WingParams wp = fit_baseline(&cf, EXT_RANGE_MIN);

    

    return 0;  
}
