#include <stdio.h>
#include "hawaii.h"

double pes(double *q) { UNUSED(q); return 0; }
void dpes(double *q, double *dVdq) { UNUSED(q); UNUSED(dVdq); }

#define N 4 

void test_fft_roundtrip()
{
    printf("---------- FFT ROUNDTRIP --------------------\n");
    double *v = (double*) malloc(N * sizeof(double));
    for (size_t i = 0; i < N; ++i) {
        v[i] = i + 1;
    }
    printf("v: ");
    for (size_t i = 0; i < N; ++i) {
       printf("%.3f ", v[i]); 
    } 
    printf("\n");

    double *fft_v = dct(v, N);
    printf("fft_v: ");
    for (size_t i = 0; i < N; ++i) {
       printf("%.3f ", fft_v[i]); 
    } 
    printf("\n");

    double *iv = idct(fft_v, N);
    
    printf("iv: ");
    for (size_t i = 0; i < N; ++i) {
       printf("%.3f ", iv[i]); 
    } 
    printf("\n");


    for (size_t i = 0; i < N; ++i) {
        assert_float_is_equal_to(iv[i], v[i], 1e-6);
    }

    free(v);
    free(fft_v);
    free(iv);
    
    printf("-------------------------------------------\n\n\n");
}

void test_pad()
{
    printf("---------- TESTING PADDING TO CLOSEST POWER OF 2 ------------\n");
    size_t len = N + 1;
    double *v = (double*) malloc(len * sizeof(double));
    memset(v, 0, len * sizeof(double));
    v[0] = 2.0;
    v[2] = 3.0;
    printf("original array (len = %zu): ", len);
    for (size_t i = 0; i < len; ++i) {
        printf("%.3f ", v[i]);
    } 
    printf("\n");

    size_t padded_len;
    double *padded = pad_to_power_of_two(v, N + 1, &padded_len);

    printf("padded (len = %zu): ", padded_len);
    for (size_t i = 0; i < padded_len; ++i) {
        printf("%.3f ", padded[i]);
    } 
    printf("\n");

    free(v);
    free(padded);
    
    printf("-----------------------------------------------\n\n");
}

void test_dct_roundtrip()
{
    const char *filename = "He-Ar/CF-He-Ar-F-300.0.txt";

    String_Builder sb = {0};
    CFnc cf = {0};

    if (!read_correlation_function(filename, &sb, &cf)) {
        printf("ERROR: could not read the file '%s'!\n", filename);
        return; 
    }

    printf("Showing head for CF:\n");
    for (size_t i = 0; i < 5; ++i) {
        printf("%.5e %.5e\n", cf.t[i], cf.data[i]);
    } 

    SFnc sf = idct_cf_to_sf(cf);

    printf("Showing head for SF:\n");
    for (size_t i = 0; i < 5; ++i) {
        printf("%.5e %.5e\n", sf.nu[i], sf.data[i]);
    } 

    CFnc cf_roundtrip = dct_sf_to_cf(sf);
    printf("Showing head for roundtrip CF:\n");
    for (size_t i = 0; i < 5; ++i) {
        printf("%.5e %.5e\n", cf_roundtrip.t[i], cf_roundtrip.data[i]);
    } 

    for (size_t i = 0; i < cf.len; ++i) {
        assert_float_is_equal_to(cf.data[i], cf_roundtrip.data[i], 1e-34);
    } 
}

int main()
{
    test_fft_roundtrip();
    test_pad();
    test_dct_roundtrip();


    return 0;
}


