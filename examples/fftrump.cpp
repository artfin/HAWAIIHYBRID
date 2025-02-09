#include "hawaii.h"
#include "loess.hpp"

#include <regex.h>

double pes(double *q) { UNUSED(q); return 0; }
void dpes(double *q, double *dVdq) { UNUSED(q); UNUSED(dVdq); }


bool write_averaged_cfs_for_CH4_CO2(double T)
{
    CFnc cf1{};
    CFnc cf2{};
    CFnc average{};
    
    String_Builder filename = {};
    bool result = true;

    {
        sb_append_cstring(&filename, "MT_CH4_CO2/p1/CF-CH4-CO2-F-");
        sb_append_format(&filename, "%.1f", T);
        sb_append_cstring(&filename, ".txt");
        if (!read_correlation_function(filename.items, NULL, &cf1)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            return_defer(false); 
        }

        sb_reset(&filename); 
    }
    { 
        sb_append_cstring(&filename, "MT_CH4_CO2/p2/CF-CH4-CO2-F-");
        sb_append_format(&filename, "%.1f", T);
        sb_append_cstring(&filename, ".txt");
        if (!read_correlation_function(filename.items, NULL, &cf2)) {
            printf("ERROR: could not read the file '%s'!\n", filename.items);
            return_defer(false); 
        }
        
        sb_reset(&filename); 
    }

    if (average_correlation_functions(&cf1, &cf2, &average) < 0) {
        return_defer(false);
    }
    
    {
        sb_append_cstring(&filename, "MT_CH4_CO2/free-final/CF-CH4-CO2-F-");
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
    free_cfnc(average);
    sb_free(&filename);

    return result; 
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

    SFnc sf = dct_numeric_sf(cf, &wp);
    sf.Temperature = T;

    double M0 = compute_M0_from_sf(sf);
    double M2 = compute_M2_from_sf(sf);
    printf("\n-------------------------------------\n");
    printf("M0 from raw SF: %.6e\n", M0);
    printf("M2 from raw SF: %.6e\n", M2);
    printf("-------------------------------------\n");

    Spectrum spraw = compute_alpha(desymmetrize_sch(sf));
    spraw.len = 25000;

    sb_reset(&filename); 
    sb_append_cstring(&filename, "MT_CH4_CO2/free-final/SPD3-CH4-CO2-F-");
    sb_append_format(&filename, "%.1f", T);
    sb_append_cstring(&filename, "-raw.txt");
    if (!writetxt(filename.items, spraw.nu, spraw.data, spraw.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename.items);
    }

    LOESS loess(sf.nu, sf.data, sf.len);
    double dnu = sf.nu[1] - sf.nu[0];

    double numin = 0.0;
    double numax = 1000.0;
    size_t npoints = 10001;
    size_t wsmin = 30;
    double wsstep = 0.66;
    size_t wsdelay = 100; 
    std::vector<double> sfsmvals = loess.eval_linear_ws(2, numin, numax, npoints, wsmin, wsstep, wsdelay, dnu); 

    free(sf.nu);
    sf.nu = linspace(numin, numax, npoints);
    memcpy(sf.data, sfsmvals.data(), npoints * sizeof(double)); 
    sf.len = npoints;

    Spectrum spd3 = compute_alpha(desymmetrize_sch(sf)); 
    
    sb_reset(&filename); 
    sb_append_cstring(&filename, "MT_CH4_CO2/free-final/SPD3-CH4-CO2-F-");
    sb_append_format(&filename, "%.1f", T);
    sb_append_cstring(&filename, ".txt");
    if (!writetxt(filename.items, spd3.nu, spd3.data, spd3.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename.items);
    }

    sb_reset(&filename); 
    sb_reset(&sb);

    free_cfnc(cf);
    // free_sfnc(sf);
    // free_spectrum(spd3);
        
    sb_free(&sb);

    return true;
}

void process_cf_for_CO_Ar()
{
    String_Builder sb{};
    CFnc cf{};

    const char *filename = "CO-Ar/CF-CO-Ar-F-300.0.txt"; 
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

    size_t EXT_RANGE_MIN = 8192; 
    WingParams wp = fit_baseline(&cf, EXT_RANGE_MIN);
    
    SFnc sf = dct_numeric_sf(cf, &wp);
    sf.Temperature = 300;

    double M0 = compute_M0_from_sf(sf);
    double M2 = compute_M2_from_sf(sf);
    printf("\n-------------------------------------\n");
    printf("M0 from raw SF: %.6e\n", M0);
    printf("M2 from raw SF: %.6e\n", M2);
    printf("-------------------------------------\n");

    Spectrum spraw = compute_alpha(desymmetrize_sch(sf));
    spraw.len = 20000;

    filename = "CO-Ar/SPD3-CO-Ar-F-300-raw.txt";
    if (!writetxt(filename, spraw.nu, spraw.data, spraw.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename);
    }

    LOESS loess(sf.nu, sf.data, sf.len);
    double dnu = sf.nu[1] - sf.nu[0];

    double numin = 0.0;
    double numax = 600.0;
    size_t npoints = 6001;
    size_t wsmin = 30;
    double wsstep = 0.7;
    size_t wsdelay = 100; 
    std::vector<double> sfsmvals = loess.eval_linear_ws(2, numin, numax, npoints, wsmin, wsstep, wsdelay, dnu); 

    free(sf.nu);
    sf.nu = linspace(numin, numax, npoints);
    memcpy(sf.data, sfsmvals.data(), npoints * sizeof(double)); 
    sf.len = npoints;

    Spectrum spd3 = compute_alpha(desymmetrize_sch(sf)); 
    
    filename = "CO-Ar/SPD3-CO-Ar-F-300.txt";
    if (!writetxt(filename, spd3.nu, spd3.data, spd3.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename);
    }

    free_cfnc(cf);
    free_sfnc(sf);
    free_spectrum(spd3);
        
    sb_free(&sb);
}


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

int main()
{
    {
        size_t ntemps = 1;
        double Temps[ntemps] = {290.0}; // , 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0};

        for (size_t i = 0; i < ntemps; ++i) {
            double T = Temps[i];
            write_averaged_cfs_for_CH4_CO2(T);
            if (!process_cfs_for_CH4_CO2(T)) {
                return 1;
            }
        }
    }
    
    // process_cf_for_CO_Ar();

    //process_correlation_function();
    // process_spectral_function();

    
    return 0;
}

