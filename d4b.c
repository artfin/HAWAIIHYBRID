#include <stdio.h>
#include "hawaii.h"

// TODO: смотрим на десимметризации при 20, 30 и 40 К, 1000 K, 2000 K, 5000 K

void desymmetrize_d4b(CFnc cf, double d0, double d1, double *m0, double *m1, Spectrum *out_spectrum, bool do_return_spectrum) 
{
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, cf.len);
    gsl_spline_init(spline, cf.t, cf.data, cf.len);

    //SFnc sf_cl = idct_cf_to_sf(cf);
    //double m0_cl = compute_Mn_from_sf_using_classical_detailed_balance(sf_cl, 0);
    //double m1_cl = compute_Mn_from_sf_using_classical_detailed_balance(sf_cl, 1);
    //printf("M0 cl = %.5e, m1 cl = %.5e\n", m0_cl, m1_cl);
    //assert(false);

    CFnc cf_d4b = {
        .t    = (double*) malloc(cf.len * sizeof(double)),
        .data = (double*) malloc(cf.len * sizeof(double)),
        .capacity = cf.len,
        .Temperature = cf.Temperature,
    };

    double cc = HBar / Boltzmann / cf.Temperature / ATU / 2.0;
    //printf("INFO: time shift constant = %.10e\n", cc);

    size_t cursor = 0;

    for (size_t i = 0; i < cf.len; ++i) {
        double t_shifted = sqrt(cf.t[i]*cf.t[i] + d1*cc*cc);

        if (t_shifted <= cf.t[cf.len - 1]) {
            cf_d4b.t[cursor] = cf.t[i];
            cf_d4b.data[cursor] = gsl_spline_eval(spline, t_shifted, acc);
            cursor++;
        }
    }
   
    cf_d4b.len = cursor;   
    //printf("INFO: len of CF after interpolation = %zu\n", cf_d4b.len);
    
    size_t padded_len = 0; 
    double *padded_t = pad_to_power_of_two(cf_d4b.t, cf_d4b.len, &padded_len); 
    double *padded_data = pad_to_power_of_two(cf_d4b.data, cf_d4b.len, &padded_len); 
    //printf("INFO: padding 't' and 'data' arrays of CF to padded_len = %zu\n", padded_len);

    free(cf_d4b.t);
    free(cf_d4b.data);
    cf_d4b.t = padded_t;
    cf_d4b.data = padded_data;
    cf_d4b.len = padded_len;

    SFnc sf_intermediate_d4b = idct_cf_to_sf(cf_d4b);
    
    //writetxt("SFD4b-F-He-Ar-50.0.txt", sf_intermediate_d4b.nu, sf_intermediate_d4b.data, sf_intermediate_d4b.len, NULL); 
    
    // for D4a
    // d0 = d0 * gsl_spline_eval(spline, 0.0, acc) / gsl_spline_eval(spline, cc, acc);

    SFnc sf_d4b = desymmetrize_schofield(sf_intermediate_d4b);
    for (size_t i = 0; i < sf_d4b.len; ++i) {
        sf_d4b.data[i] *= d0;
    }

    sf_d4b.len = 25000; 
    *m0 = compute_Mn_from_sf_using_quantum_detailed_balance(sf_d4b, 0);
    *m1 = compute_Mn_from_sf_using_quantum_detailed_balance(sf_d4b, 1);
    // printf("M0 = %.5e\n", *m0);
    // printf("M1 = %.5e\n", *m1);

    if (do_return_spectrum) {
        *out_spectrum = compute_alpha(sf_d4b); 
    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    free_cfnc(cf_d4b); 
    free_sfnc(sf_intermediate_d4b);
    free_sfnc(sf_d4b);
}

typedef struct {
    CFnc cf;
    double M0ref;
    double M1ref;
} OptimizeData;

int function_to_optimize(const gsl_vector* x, void* data, gsl_vector* f)
/*
 * [input]  x    : parameters [d0, d1]
 * [input]  data : CF
 * [output] f    : Yi - yi = difference between model and CF values
 */ 
{
    OptimizeData *opt = (OptimizeData*) data; 

    double d0 = gsl_vector_get(x, 0);
    double d1 = gsl_vector_get(x, 1);

    double m0, m1;
    desymmetrize_d4b(opt->cf, d0, d1, &m0, &m1, NULL, false);

    gsl_vector_set(f, 0, (opt->M0ref - m0) / opt->M0ref);
    gsl_vector_set(f, 1, (opt->M1ref - m1) / opt->M1ref);

    return GSL_SUCCESS;
}

void gsl_multifit_callback2(const size_t iter, void* params, const gsl_multifit_nlinear_workspace* w)
{
    UNUSED(params);

    gsl_vector* f = gsl_multifit_nlinear_residual(w);
    gsl_vector* x = gsl_multifit_nlinear_position(w);

    fprintf(stdout, "    [callback] LM iter %2zu: d0 = %.12f, d1 = %.12f, |f(x)| = %.4f\n",
            iter, gsl_vector_get(x, 0), gsl_vector_get(x, 1), gsl_blas_dnrm2(f));
}

void gsl_multifit_callback3(const size_t iter, void* params, const gsl_multifit_nlinear_workspace* w)
{
    UNUSED(params);

    gsl_vector* f = gsl_multifit_nlinear_residual(w);
    gsl_vector* x = gsl_multifit_nlinear_position(w);

    fprintf(stdout, "    [callback] LM iter %2zu: d00 = %.12f, d01 = %.12f, d02 = %.12f, d10 = %.12f, d11 = %.12f, d12 = %.12f, |f(x)| = %.4f\n",
            iter, gsl_vector_get(x, 0), gsl_vector_get(x, 1), gsl_vector_get(x, 2), 
                  gsl_vector_get(x, 3), gsl_vector_get(x, 4), gsl_vector_get(x, 5), gsl_blas_dnrm2(f));
}

void optimize_d0_and_d1(CFnc cf, double m0ref, double m1ref, double *out_d0, double *out_d1)
{
    OptimizeData opt = {
        .cf = cf,
        .M0ref = m0ref, 
        .M1ref = m1ref, 
    };

    const gsl_multifit_nlinear_type* T = gsl_multifit_nlinear_trust;

    gsl_multifit_nlinear_fdf fdf;
    fdf.f      = function_to_optimize;
    fdf.df     = NULL; // set to NULL for finite-difference Jacobian
    fdf.fvv    = NULL; // not using geodesic acceleration
    fdf.n      = 2;
    fdf.p      = 2;
    fdf.params = (void*) &opt;
    
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
    gsl_multifit_nlinear_workspace* w = gsl_multifit_nlinear_alloc(T, &fdf_params, 2, 2);
    
    double weights[] = {1.0, 1.0}; 
    gsl_vector_view wts = gsl_vector_view_array(weights, 2);

    double initial[] = {1.0, 0.1};
    gsl_vector_view initial_view = gsl_vector_view_array(initial, 2);

    gsl_multifit_nlinear_winit(&initial_view.vector, &wts.vector, &fdf, w); // initialize solver with starting point and weights

    // compute initial cost function
    gsl_vector* f = gsl_multifit_nlinear_residual(w);

    double chisq0;
    gsl_blas_ddot(f, f, &chisq0);

    int status, info;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 1e-8;

    const size_t niter_max = 100;
    status = gsl_multifit_nlinear_driver(niter_max, xtol, gtol, ftol, gsl_multifit_callback2, NULL, &info, w);

    // covariance of best fit parameters
    // gsl_matrix* J = gsl_multifit_nlinear_jac(w);
    // gsl_matrix* covar = gsl_matrix_alloc(3, 3);
    // gsl_multifit_nlinear_covar(J, 0.0, covar);

    double chisq;
    gsl_blas_ddot(f, f, &chisq); // final cost

    fprintf(stderr, "    summary from method '%s/%s'\n", gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "    number of iterations: %zu\n", gsl_multifit_nlinear_niter(w));
    fprintf(stderr, "    function evaluations: %zu\n", fdf.nevalf);
    fprintf(stderr, "    Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "    reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "    initial |f(x)| = %f\n", sqrt(chisq0));
    fprintf(stderr, "    final   |f(x)| = %f\n", sqrt(chisq));

    // double dof = n - nparams;
    // double c = fmax(1.0, sqrt(chisq / dof));

    // fprintf(stderr, "chisq / dof = %g\n", chisq / dof);
    fprintf(stderr, "d0 = %.10f\n", gsl_vector_get(w->x, 0));
    fprintf(stderr, "d1 = %.10f\n", gsl_vector_get(w->x, 1));

    *out_d0 = gsl_vector_get(w->x, 0);
    *out_d1 = gsl_vector_get(w->x, 1);
    
    gsl_multifit_nlinear_free(w);
}

void optimize_one_by_one()
{
    FILE *fp = fopen("dvals.txt", "w");
    printf("INFO: opened dvals.txt to write d0 and d1 values\n");

    //double Temperatures[] = {50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 
    //                         170.0, 190.0, 210.0, 230.0, 250.0, 270.0, 290.0};
        
    double Temperatures[] = { 1000.0 };
    size_t nTemperatures = sizeof(Temperatures)/sizeof(Temperatures[0]);

    String_Builder filename = {0};

    for (size_t i = 0; i < nTemperatures; ++i) 
    {
        printf("\n\n");

        double Temperature = Temperatures[i]; 
        CFnc cf_f = {0};
        CFnc cf_b = {0}; 

        sb_reset(&filename);
        sb_append_format(&filename, "He-Ar/CF-F/CF-F-He-Ar-%.1f.txt", Temperature); 
        if (!read_correlation_function(filename.items, NULL, &cf_f)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            exit(1); 
        }

        CFnc cf = {
            .t = cf_f.t,
            .data = malloc(cf_f.capacity * sizeof(double)),
            .len = cf_f.len,
            .capacity = cf_f.capacity,
            .Temperature = cf_f.Temperature, 
        };

        for (size_t i = 0; i < cf_f.len; ++i) {
            cf.data[i] = cf_f.data[i];
        }


        sb_reset(&filename); 
        sb_append_format(&filename, "He-Ar/CF-B/CF-B-He-Ar-%.1f.txt", Temperature); 
        if (!read_correlation_function(filename.items, NULL, &cf_b)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            // exit(1); 
        } else {
            for (size_t i = 0; i < cf_f.len; ++i) {
                cf.data[i] += cf_b.data[i];
            }
        }

        double m0ref[] = {
            //1.588844492e-06, // 50  
            //1.814755710e-06, // 70
            //2.084497471e-06, // 90
            //2.375063800e-06, // 110
            //2.677956101e-06, // 130
            //2.989174251e-06, // 150
            //3.306514410e-06, // 170
            //3.628619784e-06, // 190
            //3.954583724e-06, // 210
            //4.283761525e-06, // 230
            //4.615672169e-06, // 250
            //4.949942964e-06, // 270
            //5.286276349e-06, // 290
            //5.624428938e-06, // 310
            //5.964197698e-06, // 330
            //6.305410489e-06, // 350
            //6.647919378e-06, // 370
            //6.991595780e-06, // 390
            //7.336326838e-06, // 410
            //7.682012676e-06, // 430
            //8.028564262e-06, // 450
            //8.375901722e-06, // 470
            //8.723952993e-06, // 490, quantum
            //1.064821641e-05, // 600, hbar^4
            //1.240792216e-05, // 700, hbar^4
            //1.417303674e-05, // 800, hbar^4
            //1.594058518e-05, // 900, hbar^4
            1.770851147e-05, // 1000, hbar^4
        }; 
            
        double m1ref[] = {
            //7.537432203e-05, // 50
            //8.385039936e-05, // 70
            //9.417750094e-05, // 90
            //1.052460748e-04, // 110
            //1.166659866e-04, // 130
            //1.282654866e-04, // 150
            //1.399574642e-04, // 170
            //1.516932167e-04, // 190
            //1.634435213e-04, // 210
            //1.751898926e-04, // 230
            //1.869201569e-04, // 250
            //1.986260475e-04, // 270
            //2.103018195e-04, // 290
            //2.219434135e-04, // 310
            //2.335479297e-04, // 330
            //2.451132861e-04, // 350
            //2.566379899e-04, // 370
            //2.681209797e-04, // 390
            //2.795615160e-04, // 410
            //2.909591012e-04, // 430
            //3.023134226e-04, // 450
            //3.136243093e-04, // 470
            //3.248917000e-04, // 490
            //3.860833542e-04, // 600, hbar^4 
            //4.406031116e-04, // 700, hbar^4
            //4.940992375e-04, // 800, hbar^4
            //5.466124826e-04, // 900, hbar^4
            5.981823007e-04, // 1000, hbar^4
        };

        /*
           double d0 = 1.0;
           double d1 = 1.0;
           double m0, m1;
           Spectrum sp = desymmetrize_d4b(cf, d0, d1, &m0, &m1);
           writetxt("d4b/SPD4b-T-He-Ar-300.0.txt", sp.nu, sp.data, sp.len, NULL); 
           */

        double opt_d0, opt_d1; 
        optimize_d0_and_d1(cf, m0ref[i], m1ref[i], &opt_d0, &opt_d1);

        Spectrum sp;

        double m0, m1;
        desymmetrize_d4b(cf, opt_d0, opt_d1, &m0, &m1, &sp, true);

        sb_reset(&filename);
        sb_append_format(&filename, "d4b/SPD4b-T-He-Ar-%.1f.txt", Temperature);
        writetxt(filename.items, sp.nu, sp.data, sp.len, NULL); 

        printf("M0ref = %.8e => M0 = %.8e\n", m0ref[i], m0);
        printf("M1ref = %.8e => M1 = %.8e\n", m1ref[i], m1);
        fprintf(fp, "%.2f %.5e %.5e\n", Temperature, opt_d0, opt_d1);

        free_spectrum(sp);
        free_cfnc(cf_f);
        free_cfnc(cf_b);
    }
        
    sb_free(&filename);
    fclose(fp);
}

typedef struct {
    CFnc *cf;
    double *Temperatures;
    size_t nTemperatures;
    double *M0ref;
    double *M1ref;
} MultiTempOptimizeData;

int multitemp_optimize(const gsl_vector* x, void* data, gsl_vector* f)
/*
 * [input]  x    : parameters [d00, d01, d02, d10, d11, d12]
 * [input]  data : CF
 * [output] f    : Yi - yi = difference between model and CF values
 */ 
{
    MultiTempOptimizeData *opt = (MultiTempOptimizeData*) data; 

    double d00 = gsl_vector_get(x, 0);
    double d01 = gsl_vector_get(x, 1);
    double d02 = gsl_vector_get(x, 2);
    double d10 = gsl_vector_get(x, 3);
    double d11 = gsl_vector_get(x, 4);
    double d12 = gsl_vector_get(x, 5);

    for (size_t i = 0; i < opt->nTemperatures; ++i) {
        double T = opt->Temperatures[i];
        double d0 = d00 + d01/T + d02/T/T;
        double d1 = d10 + d11/T + d12/T/T;

        assert(T == opt->cf[i].Temperature);
        double m0, m1;
        desymmetrize_d4b(opt->cf[i], d0, d1, &m0, &m1, NULL, false);

        gsl_vector_set(f, 2*i, (opt->M0ref[i] - m0) / opt->M0ref[i]);
        gsl_vector_set(f, 2*i + 1, (opt->M1ref[i] - m1) / opt->M1ref[i]);
    }

    return GSL_SUCCESS;
}


void optimize_all_temperatures()
{
    double Temperatures[] = {50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 
                             170.0, 190.0, 210.0, 230.0, 250.0, 270.0, 290.0,
                             410.0, 430.0, 450.0, 470.0, 490.0, 
                             600.0, 700.0, 800.0, 900.0, 1000.0};

    int nTemperatures = sizeof(Temperatures)/sizeof(Temperatures[0]); 
    CFnc cf[nTemperatures];

    String_Builder filename = {};

    for (size_t i = 0; i < sizeof(Temperatures)/sizeof(Temperatures[0]); ++i) 
    {
        printf("\n\n");

        double Temperature = Temperatures[i]; 
        CFnc cf_f = {0};
        CFnc cf_b = {0}; 

        sb_reset(&filename);
        sb_append_format(&filename, "He-Ar/CF-F/CF-F-He-Ar-%.1f.txt", Temperature); 
        if (!read_correlation_function(filename.items, NULL, &cf_f)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            exit(1); 
        }

        cf[i] = (CFnc) {
            .t = malloc(cf_f.capacity * sizeof(double)), 
            .data = malloc(cf_f.capacity * sizeof(double)),
            .len = cf_f.len,
            .capacity = cf_f.capacity,
            .Temperature = cf_f.Temperature, 
        };

        for (size_t j = 0; j < cf_f.len; ++j) {
            cf[i].t[j] = cf_f.t[j];
            cf[i].data[j] = cf_f.data[j];
        }


        sb_reset(&filename); 
        sb_append_format(&filename, "He-Ar/CF-B/CF-B-He-Ar-%.1f.txt", Temperature); 
        if (!read_correlation_function(filename.items, NULL, &cf_b)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            // exit(1); 
        } else {
            for (size_t j = 0; j < cf_f.len; ++j) {
                cf[i].data[j] += cf_b.data[j];
            }
        }

        free_cfnc(cf_f);
        free_cfnc(cf_b);
    }

    double m0ref[] = {
        1.588844492e-06, // 50  
        1.814755710e-06, // 70
        2.084497471e-06, // 90
        2.375063800e-06, // 110
        2.677956101e-06, // 130
        2.989174251e-06, // 150
        3.306514410e-06, // 170
        3.628619784e-06, // 190
        3.954583724e-06, // 210
        4.283761525e-06, // 230
        4.615672169e-06, // 250
        4.949942964e-06, // 270
        5.286276349e-06, // 290
        //5.624428938e-06, // 310
        //5.964197698e-06, // 330
        //6.305410489e-06, // 350
        //6.647919378e-06, // 370
        //6.991595780e-06, // 390
        7.336326838e-06, // 410
        7.682012676e-06, // 430
        8.028564262e-06, // 450
        8.375901722e-06, // 470
        8.723952993e-06, // 490, quantum
        1.064821641e-05, // 600, hbar^4
        1.240792216e-05, // 700, hbar^4
        1.417303674e-05, // 800, hbar^4
        1.594058518e-05, // 900, hbar^4
        1.770851147e-05, // 1000, hbar^4
    }; 
        
    double m1ref[] = {
        7.537432203e-05, // 50
        8.385039936e-05, // 70
        9.417750094e-05, // 90
        1.052460748e-04, // 110
        1.166659866e-04, // 130
        1.282654866e-04, // 150
        1.399574642e-04, // 170
        1.516932167e-04, // 190
        1.634435213e-04, // 210
        1.751898926e-04, // 230
        1.869201569e-04, // 250
        1.986260475e-04, // 270
        2.103018195e-04, // 290
        //2.219434135e-04, // 310
        //2.335479297e-04, // 330
        //2.451132861e-04, // 350
        //2.566379899e-04, // 370
        //2.681209797e-04, // 390
        2.795615160e-04, // 410
        2.909591012e-04, // 430
        3.023134226e-04, // 450
        3.136243093e-04, // 470
        3.248917000e-04, // 490
        3.860833542e-04, // 600, hbar^4 
        4.406031116e-04, // 700, hbar^4
        4.940992375e-04, // 800, hbar^4
        5.466124826e-04, // 900, hbar^4
        5.981823007e-04, // 1000, hbar^4
    };

    MultiTempOptimizeData opt = {
        .cf            = cf,
        .Temperatures  = Temperatures,
        .nTemperatures = sizeof(Temperatures)/sizeof(Temperatures[0]), 
        .M0ref         = m0ref,
        .M1ref         = m1ref,
    };

    const gsl_multifit_nlinear_type* T = gsl_multifit_nlinear_trust;

    gsl_multifit_nlinear_fdf fdf;
    fdf.f      = multitemp_optimize;
    fdf.df     = NULL; // set to NULL for finite-difference Jacobian
    fdf.fvv    = NULL; // not using geodesic acceleration
    fdf.n      = 2*nTemperatures;
    fdf.p      = 6;
    fdf.params = (void*) &opt;
    
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
    gsl_multifit_nlinear_workspace* w = gsl_multifit_nlinear_alloc(T, &fdf_params, fdf.n, fdf.p);
    
    double *weights = malloc(fdf.n * sizeof(double));
    for (size_t i = 0; i < fdf.n; ++i) { 
        weights[i] = 1.0;
    }
    gsl_vector_view wts = gsl_vector_view_array(weights, fdf.n);

    double initial[] = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    gsl_vector_view initial_view = gsl_vector_view_array(initial, fdf.p);

    gsl_multifit_nlinear_winit(&initial_view.vector, &wts.vector, &fdf, w); // initialize solver with starting point and weights

    // compute initial cost function
    gsl_vector* f = gsl_multifit_nlinear_residual(w);

    double chisq0;
    gsl_blas_ddot(f, f, &chisq0);

    int status, info;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 1e-8;

    const size_t niter_max = 100;
    status = gsl_multifit_nlinear_driver(niter_max, xtol, gtol, ftol, gsl_multifit_callback3, NULL, &info, w);

    // covariance of best fit parameters
    // gsl_matrix* J = gsl_multifit_nlinear_jac(w);
    // gsl_matrix* covar = gsl_matrix_alloc(3, 3);
    // gsl_multifit_nlinear_covar(J, 0.0, covar);

    double chisq;
    gsl_blas_ddot(f, f, &chisq); // final cost

    fprintf(stderr, "    summary from method '%s/%s'\n", gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "    number of iterations: %zu\n", gsl_multifit_nlinear_niter(w));
    fprintf(stderr, "    function evaluations: %zu\n", fdf.nevalf);
    fprintf(stderr, "    Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "    reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "    initial |f(x)| = %f\n", sqrt(chisq0));
    fprintf(stderr, "    final   |f(x)| = %f\n", sqrt(chisq));

    // double dof = n - nparams;
    // double c = fmax(1.0, sqrt(chisq / dof));

    // fprintf(stderr, "chisq / dof = %g\n", chisq / dof);
    fprintf(stderr, "d00 = %.10f\n", gsl_vector_get(w->x, 0));
    fprintf(stderr, "d01 = %.10f\n", gsl_vector_get(w->x, 1));
    fprintf(stderr, "d02 = %.10f\n", gsl_vector_get(w->x, 2));
    fprintf(stderr, "d10 = %.10f\n", gsl_vector_get(w->x, 3));
    fprintf(stderr, "d11 = %.10f\n", gsl_vector_get(w->x, 4));
    fprintf(stderr, "d12 = %.10f\n", gsl_vector_get(w->x, 5));

    FILE *fp = fopen("d4b/dvals-multifit-inverseT.txt", "w");

    for (size_t i = 0; i < nTemperatures; ++i) 
    {
        double T = Temperatures[i];
        
        double d00 = gsl_vector_get(w->x, 0);
        double d01 = gsl_vector_get(w->x, 1);
        double d02 = gsl_vector_get(w->x, 2);
        double d10 = gsl_vector_get(w->x, 3);
        double d11 = gsl_vector_get(w->x, 4);
        double d12 = gsl_vector_get(w->x, 5);

        double d0 = d00 + d01/T + d02/T/T;
        double d1 = d10 + d11/T + d12/T/T;
       
        Spectrum sp = {0};

        double m0, m1;
        desymmetrize_d4b(cf[i], d0, d1, &m0, &m1, &sp, true);
        
        sb_reset(&filename);
        sb_append_format(&filename, "d4b/SPD4b-T-He-Ar-%.1f-multifit-inverseT.txt", T);
        writetxt(filename.items, sp.nu, sp.data, sp.len, NULL); 

        printf("T = %.2f :: d0 = %.5f, d1 = %.5f, M0 (ref) = %.8e, M0 = %.8e, M1 (ref) = %.8e, M1 = %.8e\n", T, d0, d1, m0ref[i], m0, m1ref[i], m1);
        fprintf(fp, "%.2f %.5f %.5f\n", T, d0, d1);

        free_spectrum(sp);
    }

    fclose(fp);
  
    for (size_t i = 0; i < nTemperatures; ++i) {
        free_cfnc(cf[i]);
    }

    sb_free(&filename);
    // *out_d0 = gsl_vector_get(w->x, 0);
    // *out_d1 = gsl_vector_get(w->x, 1);
} 

int main()
{
    optimize_one_by_one();
    // optimize_all_temperatures();

    return 0;
}
