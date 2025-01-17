#include "hawaii.h"

double pes(double *q) { UNUSED(q); return 0; }
void dpes(double *q, double *dVdq) { UNUSED(q); UNUSED(dVdq); }
    
int main()
{
    String_Builder sb = {0};
    CFnc cf = {0};

    const char *filename = "examples/CF-F-300.0.txt"; 
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
    printf("t range: %.3lf...%.3lf\n", cf.t[0], cf.t[cf.len - 1]);
    
    size_t EXT_RANGE_MIN = 8192; // far-infrared: 70K-300K
    WingParams wp = fit_baseline(&cf, EXT_RANGE_MIN);

    SFnc sf = dct_numeric_sf(cf, &wp);
    for (size_t i = 0; i < 10; ++i) {
        printf("%.5f %.5e\n", sf.nu[i], sf.data[i]);
    }
    
    double T = 300.0;
    SFnc sfd3 = desymmetrize_sch(sf, T);
    Spectrum spd3 = compute_alpha(sfd3, T); 

    filename = "sp.txt";
    if (!writetxt(filename, spd3.nu, spd3.data, spd3.len, sb.items)) {
        printf("ERROR: could not write into the file '%s'!\n", filename);
    }

    free_cfnc(cf);
    free_sfnc(sf);
    free_sfnc(sfd3);
    free_spectrum(spd3);
    free_sb(sb);

    return 0;  
}
