#include "hawaii.h"
#include "loess.hpp"

#include <map>

// NOTE: number of OMP threads is set using:
// export OMP_NUM_THREADS=4   
// TODO: rewrite parallelization logic in pthreads explicitly

// TODO: organize command-line flags to run correlation functions for different systems

bool write_averaged_cfs_for_CH4_CO2(double T)
{
    CFnc cf1{};
    CFnc cf2{};
    CFnc cf3{};
    CFnc cf4{};
    CFnc cf5{};
    CFnc cf6{};
    CFnc cf7{};
    CFnc cf8{};
    CFnc cf9{};
    CFnc average{};
    
    String_Builder filename = {};
    bool result = true;

    {
        sb_append_cstring(&filename, "MT_CH4_CO2/b-p1/CF-CH4-CO2-B-");
        sb_append_format(&filename, "%.1f", T);
        sb_append_cstring(&filename, ".txt");
        if (!read_correlation_function(filename.items, NULL, &cf1)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            return_defer(false); 
        }

        sb_reset(&filename); 
    }
    { 
        sb_append_cstring(&filename, "MT_CH4_CO2/b-p2/CF-CH4-CO2-B-");
        sb_append_format(&filename, "%.1f", T);
        sb_append_cstring(&filename, ".txt");
        if (!read_correlation_function(filename.items, NULL, &cf2)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            return_defer(false); 
        }
        
        sb_reset(&filename); 
    }
    { 
        sb_append_cstring(&filename, "MT_CH4_CO2/b-p3/CF-CH4-CO2-B-");
        sb_append_format(&filename, "%.1f", T);
        sb_append_cstring(&filename, ".txt");
        if (!read_correlation_function(filename.items, NULL, &cf3)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            return_defer(false); 
        }
        
        sb_reset(&filename); 
    }
    { 
        sb_append_cstring(&filename, "MT_CH4_CO2/b-p4/CF-CH4-CO2-B-");
        sb_append_format(&filename, "%.1f", T);
        sb_append_cstring(&filename, ".txt");
        if (!read_correlation_function(filename.items, NULL, &cf4)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            return_defer(false); 
        }
        
        sb_reset(&filename); 
    }
    { 
        sb_append_cstring(&filename, "MT_CH4_CO2/b-p5/CF-CH4-CO2-B-");
        sb_append_format(&filename, "%.1f", T);
        sb_append_cstring(&filename, ".txt");
        if (!read_correlation_function(filename.items, NULL, &cf5)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            return_defer(false); 
        }
        
        sb_reset(&filename); 
    }
    { 
        sb_append_cstring(&filename, "MT_CH4_CO2/b-p6/CF-CH4-CO2-B-");
        sb_append_format(&filename, "%.1f", T);
        sb_append_cstring(&filename, ".txt");
        if (!read_correlation_function(filename.items, NULL, &cf6)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            return_defer(false); 
        }
        
        sb_reset(&filename); 
    }
    { 
        sb_append_cstring(&filename, "MT_CH4_CO2/b-p7/CF-CH4-CO2-B-");
        sb_append_format(&filename, "%.1f", T);
        sb_append_cstring(&filename, ".txt");
        if (!read_correlation_function(filename.items, NULL, &cf7)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            return_defer(false); 
        }
        
        sb_reset(&filename); 
    }
    { 
        sb_append_cstring(&filename, "MT_CH4_CO2/b-p8/CF-CH4-CO2-B-");
        sb_append_format(&filename, "%.1f", T);
        sb_append_cstring(&filename, ".txt");
        if (!read_correlation_function(filename.items, NULL, &cf8)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            return_defer(false); 
        }
        
        sb_reset(&filename); 
    }
    { 
        sb_append_cstring(&filename, "MT_CH4_CO2/b-p9/CF-CH4-CO2-B-");
        sb_append_format(&filename, "%.1f", T);
        sb_append_cstring(&filename, ".txt");
        if (!read_correlation_function(filename.items, NULL, &cf9)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            return_defer(false); 
        }
        
        sb_reset(&filename); 
    }

    if (average_correlation_functions(&average, cf1, cf2, cf3, cf4, cf5, cf6, cf7, cf8, cf9) < 0) {
        printf("ERROR: could not average functions\n");
        return_defer(false);
    }
    
    {
        sb_append_cstring(&filename, "MT_CH4_CO2/bound-final/CF-CH4-CO2-B-");
        sb_append_format(&filename, "%.1f", T);
        sb_append_cstring(&filename, ".txt");
        FILE *fp = fopen(filename.items, "w");

        int ret = save_correlation_function(fp, average); 
        if (ret < 0) {
            printf("ERROR: could not write averaged correlation function to %s\n", filename.items); 
        } 

        printf("INFO: wrote %d characters to %s\n", ret, filename.items); 
        fclose(fp);
    }

defer:
    free_cfnc(cf1);
    free_cfnc(cf2);
    free_cfnc(cf3);
    free_cfnc(cf4);
    free_cfnc(cf5);
    free_cfnc(cf6);
    free_cfnc(cf7);
    free_cfnc(cf8);
    free_cfnc(cf9);
    free_cfnc(average);
    sb_free(&filename);

    return result; 
}

bool process_bound_cf_for_CH4_CO2(double T)
{
    const char *filename = "MT_CH4_CO2/b-zim-p1/CF-B-CH4-CO2-300.0.txt";

    String_Builder sb{};
    CFnc cf{};

    if (!read_correlation_function(filename, &sb, &cf)) {
        printf("ERROR: could not read the file '%s'!\n", filename);
        return false; 
    }
    
    printf("INFO: Loaded correlation function from %s\n", filename);
    
    for (size_t i = 0; i < cf.len; ++i) {
        cf.t[i] = cf.t[i] * ATU; 
    }
    
    printf("Number of samples: %zu\n", cf.len);
    printf("t range: %.3e...%.3e\n", cf.t[0], cf.t[cf.len - 1]);

    size_t EXT_RANGE_MIN = 8192; 
    //WingParams wp = fit_baseline(&cf, EXT_RANGE_MIN);
    SFnc sf = idct_cf_to_sf(cf);
    
    Spectrum spraw = compute_alpha(sf);
    filename = "MT_CH4_CO2/b-zim-p1/SPD3-B-CH4-CO2-300.0-raw.txt";
    if (!writetxt(filename, spraw.nu, spraw.data, spraw.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename);
    }
    
    {
        Smoothing_Config config = {
            .degree = 3, 
            .ws_min = 30,
            .ws_step = 1.0, 
            .ws_delay = 100, 
            .ws_cap = 0,
        }; 
        
        loess_init(sf.nu, sf.data, sf.len);
        loess_weight = WEIGHT_TRICUBE; 
        //loess_debug = true;

        free(sf.nu);
        free(sf.data); 
        
        size_t npoints = 4001;    
        sf.nu = loess_create_grid(0.0, 400.0, npoints); assert(sf.nu != NULL);
        sf.data = loess_apply_smoothing(&config); assert(sf.data != NULL); 
        sf.len = npoints;
    }

    double M0 = compute_Mn_from_sf_using_classical_detailed_balance(sf, 0);
    double M2 = compute_Mn_from_sf_using_classical_detailed_balance(sf, 2);
    printf("\n-------------------------------------\n");
    printf("M0 from raw SF: %.6e\n", M0);
    printf("M2 from raw SF: %.6e\n", M2);
    printf("-------------------------------------\n");

    Spectrum sp = compute_alpha(sf);
    //spraw.len = 5000;
    // printf("INFO: cutting the raw spectrum at %.3e cm-1\n", spraw.nu[spraw.len - 1]);
   
    filename = "MT_CH4_CO2/b-zim-p1/SPD3-B-CH4-CO2-300.0.txt";
    if (!writetxt(filename, sp.nu, sp.data, sp.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename);
    }


    return true;    
}

// TODO: reset for CFnc and SFnc structs to avoid constantly allocating and deallocating memory
bool process_cfs_for_CH4_CO2(double T)
{
    String_Builder sb{};
    String_Builder filename{};
    CFnc cf{};
      
    sb_append_cstring(&filename, "MT_CH4_CO2/free-final/CF-CH4-CO2-F-");
    sb_append_format(&filename, "%.1f", T);
    sb_append_cstring(&filename, ".txt");

    if (!read_correlation_function(filename.items, &sb, &cf)) {
        printf("ERROR: could not read the file '%s'!\n", filename.items);
        return false; 
    }
    
    printf("INFO: Loaded correlation function from %s\n", filename.items);
    
    for (size_t i = 0; i < cf.len; ++i) {
        cf.t[i] = cf.t[i] * ATU; 
    }

    printf("Number of samples: %zu\n", cf.len);
    printf("t range: %.3e...%.3e\n", cf.t[0], cf.t[cf.len - 1]);

    for (size_t j = 0; j < 10; ++j) {
        printf("t = %.3e, cf value = %.5e\n", cf.t[j], cf.data[j]);
    }

    size_t EXT_RANGE_MIN = 6192; 
    WingParams wp = fit_baseline(&cf, EXT_RANGE_MIN);

    SFnc sf = idct_cf_to_sf(cf);

    double M0 = compute_Mn_from_sf_using_classical_detailed_balance(sf, 0);
    double M2 = compute_Mn_from_sf_using_classical_detailed_balance(sf, 2);
    printf("\n-------------------------------------\n");
    printf("M0 from raw SF: %.6e\n", M0);
    printf("M2 from raw SF: %.6e\n", M2);
    printf("-------------------------------------\n");

    Spectrum spraw = compute_alpha(desymmetrize_schofield(sf));
    spraw.len = 25000;
    printf("INFO: cutting the raw spectrum at %.3e cm-1\n", spraw.nu[spraw.len - 1]);

    sb_reset(&filename); 
    sb_append_cstring(&filename, "MT_CH4_CO2/free-final/SPD3-CH4-CO2-F-");
    sb_append_format(&filename, "%.1f", T);
    sb_append_cstring(&filename, "-raw.txt");
    if (!writetxt(filename.items, spraw.nu, spraw.data, spraw.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename.items);
    }

    {
        // sadly, capping the window size is not beneficial in this case. 
        // high-frequency noise is still present at ws_cap = 3500 
        // 
        // ws_step = 1.0 produces good-looking result at 200 K
        Smoothing_Config config = {
            .degree = 3, 
            .ws_min = 30,
            .ws_step = 1.0, 
            .ws_delay = 100, 
            .ws_cap = 0,
        }; 
        
        loess_init(sf.nu, sf.data, sf.len);
        loess_weight = WEIGHT_TRICUBE; 
        //loess_debug = true;

        free(sf.nu);
        free(sf.data); 
        
        size_t npoints = 10001;    
        sf.nu = loess_create_grid(0.0, 1000.0, npoints); assert(sf.nu != NULL);
        sf.data = loess_apply_smoothing(&config); assert(sf.data != NULL); 
        sf.len = npoints;
    }

    Spectrum spd3 = compute_alpha(desymmetrize_schofield(sf)); 
    
    sb_reset(&filename); 
    sb_append_cstring(&filename, "MT_CH4_CO2/free-final/SPD3-CH4-CO2-F-");
    sb_append_format(&filename, "%.1f", T);
    sb_append_cstring(&filename, ".txt");
    if (!writetxt(filename.items, spd3.nu, spd3.data, spd3.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename.items);
    }

    free_cfnc(cf);
    free_sfnc(sf);
    free_spectrum(spraw);
    free_spectrum(spd3);
        
    loess_free();
    sb_free(&sb);

    return true;
}

void process_cf_for_CO_Ar()
{
    String_Builder sb{};
    CFnc cf{};

    const char *filename = "CO-Ar/CF-CO-Ar-cross-F-300.0.txt"; 
    if (!read_correlation_function(filename, &sb, &cf)) {
        printf("ERROR: could not read the file '%s'!\n", filename);
        exit(1);
    }
    
    printf("INFO: Loaded correlation function from %s\n", filename);
    
    for (size_t i = 0; i < cf.len; ++i) {
        cf.t[i] = cf.t[i] * ATU; 
    }

    printf("Number of samples: %zu\n", cf.len);
    printf("t range: %.3e...%.3e\n", cf.t[0], cf.t[cf.len - 1]);

    for (size_t j = 0; j < 10; ++j) {
        printf("t = %.3e, cf value = %.5e\n", cf.t[j], cf.data[j]);
    }

    // size_t EXT_RANGE_MIN = 4192; 
    // WingParams wp = fit_baseline(&cf, EXT_RANGE_MIN);

    SFnc sf = idct_cf_to_sf(cf);
    
    filename = "CO-Ar/SF-CO-Ar-cross-F-300.0.txt.2"; 
    if (!writetxt(filename, sf.nu, sf.data, sf.len, NULL)) {
        printf("ERROR: could not write into the file '%s'!\n", filename);
    }

    double M0 = compute_Mn_from_sf_using_classical_detailed_balance(sf, 0);
    double M2 = compute_Mn_from_sf_using_classical_detailed_balance(sf, 2);
    printf("\n-------------------------------------\n");
    printf("M0 from raw SF: %.6e\n", M0);
    printf("M2 from raw SF: %.6e\n", M2);
    printf("-------------------------------------\n");

    Spectrum spraw = compute_alpha(sf);
    spraw.len = 20000;

    filename = "CO-Ar/SP-CO-Ar-cross-F-300.0-raw.txt.2";
    if (!writetxt(filename, spraw.nu, spraw.data, spraw.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename);
    }
   
    /* 
    {
        Smoothing_Config config = {
            .degree = 3, 
            .ws_min = 30,
            .ws_step = 0.7, 
            .ws_delay = 100, 
            .ws_cap = 0,
        }; 
        
        loess_init(sf.nu, sf.data, sf.len);
        loess_weight = WEIGHT_TRICUBE; 
        //loess_debug = true;

        free(sf.nu);
        free(sf.data); 
        
        size_t npoints = 6001;    
        sf.nu = loess_create_grid(0.0, 600.0, npoints); assert(sf.nu != NULL);
        sf.data = loess_apply_smoothing(&config); assert(sf.data != NULL); 
        sf.len = npoints;
    }

    Spectrum spd3 = compute_alpha(desymmetrize_sch(sf)); 
    
    filename = "CO-Ar/SPD3-CO-Ar-F-152.0.txt";
    if (!writetxt(filename, spd3.nu, spd3.data, spd3.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename);
    }
    free_cfnc(cf);
    free_sfnc(sf);
    free_spectrum(spd3);
        
    sb_free(&sb);
    */
}

bool write_averaged_cfs_for_N2_Ar(double T)
{
    CFnc cf1{};
    CFnc cf2{};
    CFnc cf3{};
    CFnc average{};

    String_Builder filename = {};
    bool result = true;
    
    //{
    //    sb_append_format(&filename, "MT_N2_Ar/f-p3/CF-N2-Ar-F-%.1f.txt", T);
    //    if (!read_correlation_function(filename.items, NULL, &cf1)) {
    //        printf("ERROR: could not read the file '%s'!\n", filename.items);
    //        return_defer(false); 
    //    }

    //    sb_reset(&filename); 
    //}

    //{
    //    sb_append_format(&filename, "MT_N2_Ar/f-p4/CF-N2-Ar-F-%.1f.txt", T);
    //    if (!read_correlation_function(filename.items, NULL, &cf2)) {
    //        printf("ERROR: could not read the file '%s'!\n", filename.items);
    //        return_defer(false); 
    //    }

    //    sb_reset(&filename); 
    //}

    {
        sb_append_format(&filename, "MT_N2_Ar/b-p2/CF-N2-Ar-B-%.1f.txt", T);
        if (!read_correlation_function(filename.items, NULL, &cf3)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            return_defer(false); 
        }

        sb_reset(&filename); 
    }
    
    if (average_correlation_functions(&average, cf3) < 0) {
        printf("ERROR: could not average functions\n");
        return_defer(false);
    }
    
    {
        sb_append_format(&filename, "MT_N2_Ar/bound-final/CF-N2-Ar-B-%.1f.txt");
        FILE *fp = fopen(filename.items, "w");

        int ret = save_correlation_function(fp, average); 
        if (ret < 0) {
            printf("ERROR: could not write averaged correlation function to %s\n", filename.items); 
        } 

        printf("INFO: wrote %d characters to %s\n", ret, filename.items); 
        fclose(fp);
    }

defer:
    free_cfnc(cf1);
    free_cfnc(cf2);
    free_cfnc(average);
    sb_free(&filename);

    return result;
}

void process_bound_cf_for_N2_Ar(double T)
{
    String_Builder sb{};
    String_Builder filename{};
    CFnc cf{};
    
    sb_append_format(&filename, "MT_N2_Ar/bound-final/CF-N2-Ar-B-%.1f.txt", T);
    
    if (!read_correlation_function(filename.items, &sb, &cf)) {
        printf("ERROR: could not read the file '%s'!\n", filename.items);
        exit(1);
    }
    
    printf("INFO: Loaded correlation function from %s\n", filename.items);
    
    for (size_t i = 0; i < cf.len; ++i) {
        cf.t[i] = cf.t[i] * ATU; 
    }
    
    printf("Number of samples: %zu\n", cf.len);
    printf("t range: %.3e...%.3e\n", cf.t[0], cf.t[cf.len - 1]);

    for (size_t j = 0; j < 10; ++j) {
        printf("t = %.3e, cf value = %.5e\n", cf.t[j], cf.data[j]);
    }
   
    size_t EXT_RANGE_MIN = 8192; 
    // WingParams wp = fit_baseline(&cf, EXT_RANGE_MIN);
    SFnc sf = idct_cf_to_sf(cf);
  
    sb_reset(&filename);
    sb_append_format(&filename, "MT_N2_Ar/bound-final/SF-N2-Ar-B-%.1f.txt", T);

    if (!writetxt(filename.items, sf.nu, sf.data, sf.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename.items);
    }
    
    Spectrum spraw = compute_alpha(sf);
    
    sb_reset(&filename);
    sb_append_format(&filename, "MT_N2_Ar/bound-final/SPCL-N2-Ar-B-%.1f-raw.txt", T);
    if (!writetxt(filename.items, spraw.nu, spraw.data, spraw.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename.items);
    }
    
    {
        Smoothing_Config config = {
            .degree   = 3,
            .ws_min   = 30,
            .ws_step  = 0.3,
            .ws_delay = 100,
            .ws_cap   = 0,
        }; 
        
        loess_init(sf.nu, sf.data, sf.len);
        loess_weight = WEIGHT_TRICUBE; 
        //loess_debug = true;

        free(sf.nu);
        free(sf.data); 
        
        size_t npoints = 4001;    
        sf.nu = loess_create_grid(0.0, 400.0, npoints); assert(sf.nu != NULL);
        sf.data = loess_apply_smoothing(&config);       assert(sf.data != NULL); 
        sf.len = npoints;
    }
        
    double M0 = compute_Mn_from_sf_using_classical_detailed_balance(sf, 0);
    double M2 = compute_Mn_from_sf_using_classical_detailed_balance(sf, 2);
        
    printf("\n-------------------------------------\n");
    printf("M0 from classical SF: %.6e\n", M0);
    printf("M2 from classical SF: %.6e\n", M2);
    printf("-------------------------------------\n");
       
    Spectrum spcl = compute_alpha(sf);

    sb_reset(&filename);
    sb_append_format(&filename, "MT_N2_Ar/bound-final/SPCL-N2-Ar-B-%.1f.txt", T);

    if (!writetxt(filename.items, spcl.nu, spcl.data, spcl.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename.items);
    }
    
    Spectrum spd3 = compute_alpha(desymmetrize_schofield(sf)); 
    
    sb_reset(&filename);
    sb_append_format(&filename, "MT_N2_Ar/bound-final/SPD3-N2-Ar-B-%.1f.txt", T);
    
    if (!writetxt(filename.items, spd3.nu, spd3.data, spd3.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename.items);
    }

    free_cfnc(cf);
    free_sfnc(sf);
    free_spectrum(spcl);
    free_spectrum(spd3);
        
    sb_free(&sb);
}

void process_cf_for_N2_Ar(double T)
{
    String_Builder sb{};
    String_Builder filename{};
    CFnc cf{};

    sb_append_format(&filename, "MT_N2_Ar/free-final/CF-N2-Ar-F-%.1f.txt", T);

    if (!read_correlation_function(filename.items, &sb, &cf)) {
        printf("ERROR: could not read the file '%s'!\n", filename.items);
        exit(1);
    }

    printf("INFO: Loaded correlation function from %s\n", filename.items);

    std::map<double, double> pff = {
      {500.0, 1.0},
      {490.0, 1.0},
      {480.0, 1.0},
      {470.0, 1.0},
      {460.0, 1.0},
      {450.0, 1.0},
      {440.0, 1.0},
      {430.0, 1.0},
      {420.0, 1.0},
      {410.0, 1.0},
      {400.0, 2.67966e+05 / 5.50081e+05},
      {390.0, 2.67973e+05 / 5.16356e+05},
      {380.0, 2.67991e+05 / 4.83923e+05},
      {370.0, 2.68029e+05 / 4.52775e+05},
      {360.0, 2.68097e+05 / 4.22907e+05},
      {350.0, 2.68095e+05 / 3.94145e+05},
      {340.0, 2.68166e+05 / 3.66689e+05},
      {330.0, 2.68150e+05 / 3.40299e+05},
      {320.0, 2.68244e+05 / 3.15211e+05},
      {310.0, 2.68219e+05 / 2.91132e+05},
      {300.0, 1.0},
      {290.0, 2.68307e+05 / 2.46503e+05}, 
      {280.0, 2.68286e+05 / 2.25782e+05}, 
      {270.0, 2.68408e+05 / 2.06254e+05}, 
      {260.0, 2.68455e+05 / 1.87716e+05}, 
      {250.0, 2.68456e+05 / 1.70185e+05}, 
      {240.0, 2.68496e+05 / 1.53696e+05}, 
      {230.0, 2.68582e+05 / 1.38227e+05}, 
      {220.0, 2.68659e+05 / 1.23724e+05}, 
      {210.0, 2.68785e+05 / 1.10192e+05}, 
      {200.0, 1.0},
      {190.0, 1.0}, // 2.69040e+05 / 8.58810e+04},
      {180.0, 1.0}, // 2.69078e+05 / 7.50338e+04},
      {170.0, 1.0}, // 2.69266e+05 / 6.50880e+04},
      {160.0, 1.0}, // 2.69405e+05 / 5.59632e+04},
      {150.0, 1.0}, // 2.69600e+05 / 4.76590e+04},
      {140.0, 1.0}, // 2.69773e+05 / 4.01343e+04},
      {130.0, 1.0}, // 2.69958e+05 / 3.33696e+04},
      {120.0, 1.0}, // 2.70302e+05 / 2.73527e+04},
      {110.0, 1.0}, // 2.70594e+05 / 2.20291e+04},
      {100.0, 1.0}, // 2.70950e+05 / 1.73815e+04},
      { 90.0, 1.0}, // 2.71430e+05 / 1.33802e+04},
      { 80.0, 1.0}, // 2.72051e+05 / 9.99014e+03},
      { 70.0, 1.0}, // 2.72888e+05 / 7.17673e+03},
    };

    for (size_t i = 0; i < cf.len; ++i) {
        cf.t[i] = cf.t[i] * ATU; 
        cf.data[i] = cf.data[i] * pff[T];
    }
    
    printf("Number of samples: %zu\n", cf.len);
    printf("t range: %.3e...%.3e\n", cf.t[0], cf.t[cf.len - 1]);

    for (size_t j = 0; j < 10; ++j) {
        printf("t = %.3e, cf value = %.5e\n", cf.t[j], cf.data[j]);
    }
   
    size_t EXT_RANGE_MIN = 8192; 
    // WingParams wp = fit_baseline(&cf, EXT_RANGE_MIN);
    SFnc sf = idct_cf_to_sf(cf);
  
    sb_reset(&filename);
    sb_append_format(&filename, "MT_N2_Ar/free-final/SF-N2-Ar-F-%.1f.txt", T);

    if (!writetxt(filename.items, sf.nu, sf.data, sf.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename.items);
    }
    
    Spectrum spraw = compute_alpha(sf);

    sb_reset(&filename);
    sb_append_format(&filename, "MT_N2_Ar/free-final/SPCL-N2-Ar-F-%.1f-raw.txt", T);
    if (!writetxt(filename.items, spraw.nu, spraw.data, spraw.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename.items);
    }
    
    {
        Smoothing_Config config = {
            .degree   = 3,
            .ws_min   = 30,
            .ws_step  = 0.7,
            .ws_delay = 100,
            .ws_cap   = 0,
        }; 
        
        loess_init(sf.nu, sf.data, sf.len);
        loess_weight = WEIGHT_TRICUBE; 
        //loess_debug = true;

        free(sf.nu);
        free(sf.data); 
        
        size_t npoints = 6001;    
        sf.nu = loess_create_grid(0.0, 600.0, npoints); assert(sf.nu != NULL);
        sf.data = loess_apply_smoothing(&config);       assert(sf.data != NULL); 
        sf.len = npoints;
    }
        
    double M0 = compute_Mn_from_sf_using_classical_detailed_balance(sf, 0);
    double M2 = compute_Mn_from_sf_using_classical_detailed_balance(sf, 2);
        
    printf("\n-------------------------------------\n");
    printf("M0 from classical SF: %.6e\n", M0);
    printf("M2 from classical SF: %.6e\n", M2);
    printf("-------------------------------------\n");
       
    Spectrum spcl = compute_alpha(sf);

    sb_reset(&filename);
    sb_append_format(&filename, "MT_N2_Ar/free-final/SPCL-N2-Ar-F-%.1f.txt", T);

    if (!writetxt(filename.items, spcl.nu, spcl.data, spcl.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename.items);
    }
    
    Spectrum spd3 = compute_alpha(desymmetrize_schofield(sf)); 
    
    sb_reset(&filename);
    sb_append_format(&filename, "MT_N2_Ar/free-final/SPD3-N2-Ar-F-%.1f.txt", T);
    
    if (!writetxt(filename.items, spd3.nu, spd3.data, spd3.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename.items);
    }

    free_cfnc(cf);
    free_sfnc(sf);
    free_spectrum(spcl);
    free_spectrum(spd3);
        
    sb_free(&sb);
}

void process_cf_for_He_Ar(double T)
{
    String_Builder sb{};
    CFnc cf{};
    
    String_Builder filename{};
    sb_append_format(&filename, "He-Ar/CF-F/CF-F-He-Ar-%.1f.txt", T);

    if (!read_correlation_function(filename.items, &sb, &cf)) {
        printf("ERROR: could not read the file '%s'!\n", filename.items);
        exit(1);
    }
    
    printf("INFO: Loaded correlation function from %s\n", filename.items);
    
    printf("Number of samples: %zu\n", cf.len);
    printf("t range: %.3e...%.3e\n", cf.t[0], cf.t[cf.len - 1]);

    for (size_t j = 0; j < 10; ++j) {
        printf("t = %.3e, cf value = %.5e\n", cf.t[j], cf.data[j]);
    }
    
    // size_t EXT_RANGE_MIN = 2048; 
    // WingParams wp = fit_baseline(&cf, EXT_RANGE_MIN);
    
    SFnc sf = idct_cf_to_sf(cf);
   
    sb_reset(&filename);
    sb_append_format(&filename, "He-Ar/SF-F/SF-F-He-Ar-%.1f.txt", T);

    if (!writetxt(filename.items, sf.nu, sf.data, sf.len, NULL)) {
        printf("ERROR: could not write into the file '%s'!\n", filename.items);
    }
    
    {
        sf.len = 30000; 
        double M0 = compute_Mn_from_sf_using_classical_detailed_balance(sf, 0);
        double M1 = compute_Mn_from_sf_using_classical_detailed_balance(sf, 1);
        double M2 = compute_Mn_from_sf_using_classical_detailed_balance(sf, 2);

        printf("\n-------------------------------------\n");
        printf("M0 from raw classical SF: %.6e\n", M0);
        printf("M1 from raw classical SF: %.6e\n", M1);
        printf("M2 from raw classical SF: %.6e\n", M2);
        printf("-------------------------------------\n");

        Spectrum spcl = compute_alpha(sf);

        sb_reset(&filename);
        sb_append_format(&filename, "He-Ar/SF-F/SPCL-F-He-Ar-%.1f-raw.txt", T);
        if (!writetxt(filename.items, spcl.nu, spcl.data, spcl.len, sb.items)) {
            printf("ERROR: could not write into the file '%s'!\n", filename.items);
        }
    }

    {
        SFnc sfd3 = desymmetrize_schofield(sf);
        Spectrum spd3 = compute_alpha(sfd3);
        spd3.len = 30000;
        
        sb_reset(&filename);
        sb_append_format(&filename, "He-Ar/SF-F/SPD3-F-He-Ar-%.1f-raw.txt", T);
        if (!writetxt(filename.items, spd3.nu, spd3.data, spd3.len, sb.items)) {
            printf("ERROR: could not write into the file '%s'!\n", filename.items);
        }
    } 
}

/*
void process_correlation_function() 
{
    double Temperature = 270.0;

    String_Builder sb{};
    CFnc cf{};

    //const char *filename = "examples/CF-F-300.0.txt"; 
    //const char *filename = "./CH4-CO2/CF-CH4-CO2-F-300.0.txt"; 
    // const char *filename = "./CO-Ar/CF-CO-AR-300.0.txt"; 
    
    const char *filename = "MT_CO2_Ar/CF-F-270.0.txt";
    if (!read_correlation_function(filename, &sb, &cf)) {
        printf("ERROR: could not read the file '%s'!\n", filename);
        exit(1);
    } 
   
    printf("INFO: Loaded correlation function from %s\n", filename);
    printf("HEADER:\n %s\n\n", sb.items);
    
    for (size_t i = 0; i < cf.len; ++i) {
        cf.t[i] = cf.t[i] * ATU; 
    }

    printf("Number of samples: %zu\n", cf.len);
    printf("t range: %.3e...%.3e\n", cf.t[0], cf.t[cf.len - 1]);
    
    size_t EXT_RANGE_MIN = 8192; // far-infrared CO2-Ar: 70K-300K
    WingParams wp = fit_baseline(&cf, EXT_RANGE_MIN);
    // WingParams wp = { .A = 0.0, .B = 0.0, .C = 0.0 };

    SFnc sf = dct_numeric_sf(cf, &wp);

    for (size_t i = 0; i < 10; ++i) {
        printf("%.5f %.5e\n", sf.nu[i], sf.data[i]);
    }

    double M0 = compute_M0_from_sf(sf);
    double M2 = compute_M2_from_sf(sf);
    printf("\n-------------------------------------\n");
    printf("M0 from raw SF: %.6e\n", M0);
    printf("M2 from raw SF: %.6e\n", M2);
    printf("-------------------------------------\n");

    Spectrum spraw = compute_alpha(desymmetrize_sch(sf));
    spraw.len = 20000;  
    filename = "MT_CO2_Ar/sp_raw.txt";
    if (!writetxt(filename, spraw.nu, spraw.data, spraw.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename);
    }
   
    LOESS loess(sf.nu, sf.data, sf.len);
    double dnu = sf.nu[1] - sf.nu[0];

    double numin = 0.0;
    double numax = 300.0;
    size_t npoints = 3001;
    size_t wsmin = 30;
    double wsstep = 1.0;
    size_t wsdelay = 100; 
    std::vector<double> sfsmvals = loess.eval_linear_ws(2, numin, numax, npoints, wsmin, wsstep, wsdelay, dnu); 

    free(sf.nu);
    sf.nu = linspace(numin, numax, npoints);
    memcpy(sf.data, sfsmvals.data(), npoints * sizeof(double)); 
    sf.len = npoints;

    Spectrum spd3 = compute_alpha(desymmetrize_sch(sf)); 

    filename = "MT_CO2_Ar/SPD3-F-270.0.txt";
    if (!writetxt(filename, spd3.nu, spd3.data, spd3.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename);
    }

    free_cfnc(cf);
    free_sfnc(sf);
    free_spectrum(spd3);
    sb_free(&sb);
}

void process_spectral_function()
{
    double Temperature = 300.0;

    String_Builder sb{};
    memset(&sb, 0, sizeof(sb));

    SFnc sf{};

    const char *filename = "examples/SF-PRMU-300.0.txt"; 
    if (!read_spectral_function(filename, &sb, &sf)) {
        printf("ERROR: could not read the file '%s'!\n", filename);
        exit(1);
    } 
   
    printf("INFO: Loaded spectral function from %s\n", filename);
    printf("HEADER:\n %s\n\n", sb.items);
    
    printf("Number of samples: %zu\n", sf.len);
    printf("Frequency range: %.3lf...%.3lf\n", sf.nu[0], sf.nu[sf.len - 1]);

    sf.Temperature = Temperature;    
    Spectrum spd3 = compute_alpha(desymmetrize_sch(sf)); 
  
    filename = "SP-PRMU-300.0.txt";
    if (!writetxt(filename, spd3.nu, spd3.data, spd3.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'\n", filename);
    }

    free_sfnc(sf);
    free_spectrum(spd3);
    sb_free(&sb);
}
*/

int main()
{
    // He-Ar
    process_cf_for_He_Ar(1000.0);

    // N2-Ar
    // {
    //     double T = 150.0;
    //     write_averaged_cfs_for_N2_Ar(T);
    //     process_bound_cf_for_N2_Ar(T);
    //     //process_cf_for_N2_Ar(T);
    // }

    // CH4-CO2
    //write_averaged_cfs_for_CH4_CO2(300.0);
    //process_bound_cf_for_CH4_CO2(300.0);
    /*
    {
        size_t ntemps = 1;
        double Temps[ntemps] = {230.0};

        for (size_t i = 0; i < ntemps; ++i) {
            double T = Temps[i];
            write_averaged_cfs_for_CH4_CO2(T);
            if (!process_cfs_for_CH4_CO2(T)) {
                return 1;
            }
        }
    }
    */

    // process_cf_for_CO_Ar();

    //process_correlation_function();
    // process_spectral_function();

    
    return 0;
}

