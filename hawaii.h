#ifndef HAWAII_H_
#define HAWAII_H_

#include <assert.h>
#include <errno.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <regex.h>
#include <time.h>
#include <unistd.h>

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_spline.h>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifndef __cplusplus
#define _GNU_SOURCE
void sincos(double, double*, double*);
int syncfs(int);
#endif
#include <math.h>
#include <complex.h>

#ifndef __cplusplus
#define MT_GENERATE_CODE_IN_HEADER 0
#endif 
#include "mtwist.h"
/*
 * We are using the following functions from mtwist:
 * double mt_drand(void)
 *   Return a pseudorandom double in [0,1) with 32 bits of randomness
 *
 * uint32_t mt_lrand(void);
 *   Generate 32-bit random value 
 */


#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>

#include "array.h"
#include "constants.h"

#include "arena.h"

#define IPHI    0
#define IPPHI   1
#define ITHETA  2
#define IPTHETA 3
#define IR      4
#define IPR     5
// used only in Monomer.qp
#define IPSI    4
#define IPPSI   5

#define UNUSED(x) (void)(x)
#define TODO(message) do { fprintf(stderr, "%s:%d: TODO: %s\n", __FILE__, __LINE__, message); abort(); } while(0)
#define UNREACHABLE(message) do { fprintf(stderr, "%s:%d: UNREACHABLE: %s\n", __FILE__, __LINE__, message); abort(); } while(0)

#define return_defer(value) do { result = (value); goto defer; } while(0)

#define DA_INIT_CAP 256
#define da_append(da, item)                                                          \
    do {                                                                             \
        if ((da)->count >= (da)->capacity) {                                         \
            (da)->capacity = (da)->capacity == 0 ? DA_INIT_CAP : (da)->capacity*2;   \
            (da)->items = realloc((da)->items, (da)->capacity*sizeof(*(da)->items)); \
            assert((da)->items != NULL && "ASSERT: not enough memory\n");            \
        }                                                                            \
                                                                                     \
        (da)->items[(da)->count++] = (item);                                         \
    } while (0)

#define da_insert(da, i, item) \
    do {                                 \
        if ((i < 0) || ((i) > (da)->count)) { \
            assert(0 && "ASSERT: index out of bounds\n"); \
        } \
        if ((da)->count >= (da)->capacity) { \
            (da)->capacity = (da)->capacity == 0 ? DA_INIT_CAP : (da)->capacity*2; \
            (da)->items = realloc((da)->items, (da)->capacity*sizeof(*(da)->items)); \
            assert((da)->items != NULL && "ASSERT: not enough memory\n"); \
        } \
        memmove((da)->items + (i) + 1, (da)->items + (i), ((da)->count - (i))*sizeof(*(da)->items)); \
        (da)->items[(i)] = (item); \
        (da)->count++; \
    } while(0)

#ifdef USE_MPI
#define INIT_WRANK                          \
    int _wrank;                             \
    MPI_Comm_rank(MPI_COMM_WORLD, &_wrank); \

#define INIT_WSIZE                          \
    int _wsize;                             \
    MPI_Comm_size(MPI_COMM_WORLD, &_wsize); \

#define PRINT0(...) if (_wrank == 0) printf(__VA_ARGS__); 
#else
#define INIT_WRANK
#define INIT_WSIZE
#define PRINT0(...) printf(__VA_ARGS__)
#endif

#ifdef __cplusplus
extern "C" {
#endif

// This enum allows us to both differentiate between the systems of different type
// as well as to store the size of phase point: 
//   size(phase_point) = MonomerType%MODULO_BASE
#define MODULO_BASE 100
typedef enum {
    ATOM                                 = 0,
    LINEAR_MOLECULE                      = 4,
    LINEAR_MOLECULE_REQ_INTEGER          = MODULO_BASE + 4,
    LINEAR_MOLECULE_REQ_HALFINTEGER      = 2*MODULO_BASE + 4,
    ROTOR                                = 6,
    ROTOR_REQUANTIZED_ROTATION           = MODULO_BASE + 6,
} MonomerType;

// 17.01.2025 NOTE: changed 'double I[3]' -> 'double II[3]' to avoid collision when including <complex.h> 
typedef struct { 
    MonomerType t;
    double II[3];
    double *qp;
    double *dVdq;
    // TODO: should this be a monomer field or a bool[2] in MoleculeSystem?
    bool apply_requantization;
} Monomer; 

typedef struct {
    double intermolecular_qp[6];
    Monomer m1;
    Monomer m2;
    double mu;

    size_t Q_SIZE;
    size_t QP_SIZE;

    double *intermediate_q;
    double *dVdq;

    size_t seed; 

    // time_t: time as the number of seconds since the Epoch (1970-01-01)  
    time_t init_rawtime; // set in init_ms 
    time_t temp_rawtime; // used to mark the previous iteration of a longlasting calculation
} MoleculeSystem;

typedef enum {
    FREE_AND_METASTABLE, // автоматически отделяем вклады
    BOUND,
} PairState;

typedef struct {
    PairState ps;

    /* sampling */
    double sampler_Rmin; // a.u.
    double sampler_Rmax; // a.u.
    double pesmin; // Hartree 

    /* initial spectral moments check */ 
    size_t initialM0_npoints;
    size_t initialM2_npoints;
    double partial_partition_function_ratio;
   
    /* requantization */ 
    size_t torque_cache_len;
    double torque_limit;

    /* trajectory */
    double sampling_time; // a.t.u.
    size_t MaxTrajectoryLength;
    double cvode_tolerance;
   
    /* applicable to both correlation function AND spectral function calculations */ 
    size_t niterations; // the accumulation of the total_trajectories is split into this number of iterations
    size_t total_trajectories; // the total number of accumulated trajectories (ALL iterations) 
    
    /* correlation calculation ONLY */
    const char* cf_filename;
    double Rcut; // distance at which the trajectory is forcefully stopped, a.u. 

    /* pr/mu calculation ONLY */   
    const char *sf_filename;
    double ApproximateFrequencyMax; // cm-1
    double R0; // initial distance, a.u. 
    double average_time_between_collisions; // a.t.u.

    /* correlation function array ONLY */
    double* partial_partition_function_ratios;
    double *satellite_temperatures;
    size_t num_satellite_temperatures;
    const char **cf_filenames;
} CalcParams;


typedef enum {
    PRMU,
    CORRELATION_SINGLE,
    CORRELATION_ARRAY,
} CalculationType;

typedef enum {
    REJECT,
    MCMC
} SamplingType;


typedef struct {
    double *t; 
    double *data;
    size_t len;      // # of samples in *t, *data
    size_t capacity; // capacity *t, *data
    double ntraj;    // # of trajectories used for averaging
    double Temperature;
    bool normalized; // flag that indicates whether the *data samples are normalized by # of trajectories 
} CFnc;

/*
 * An array of correlation functions for a fixed 'base temperature'
 * 'base temperature'      -- a temperature for which sampling of initial conditions is performed
 * 'satellite temperature' -- a target temperature for which correlation function is re-weighted 
 */ 
typedef struct {
    double *t;
    double **data;
    size_t ntemp;
    size_t len;    // # of samples in *t and elements of *data 
    double *nstar; // array of # effective trajectories
    size_t ntraj; // # of calculated trajectories
} CFncArray;

typedef struct {
    double *nu;
    double *data;
    size_t len;      // # of samples in *nu, *data
    size_t capacity; // capacity of *nu, *data
    size_t ntraj;    // # of trajectories used for averaging
    double Temperature;       
} SFnc;

typedef struct {
    double *nu;
    double *data;
    size_t len;      // # of samples of *nu, *data
    size_t capacity; // capacity of *nu, *data
} Spectrum;

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

typedef struct {
    double before2;
    double before;
    double current;

    size_t turning_points;
    
    size_t called;
    bool ready;
} Tracker;


/* ----------------------------- */
/*    User-Supplied Functions    */
/* ------------------------------*/
typedef void (*dipolePtr)(double*, double[3]);
extern dipolePtr dipole;

extern double pes(double *q);
extern void dpes(double *q, double *dVdq);
/* ------------------------------*/

MoleculeSystem init_ms(double mu, MonomerType t1, MonomerType t2, double *I1, double *I2, size_t seed); 
void free_ms(MoleculeSystem *ms);

const char* monomer_type_name(MonomerType t);
const char* pair_state_name(PairState ps);

void put_qp_into_ms(MoleculeSystem *ms, Array qp);
void get_qp_from_ms(MoleculeSystem *ms, Array *qp);

// this function takes the phase point from MoleculeSystem 
// and packs it into "intermediate_q" field. Then this pointer
// can be passed to potential energy / its derivatives / dipole function 
void extract_q_and_write_into_ms(MoleculeSystem *ms);

// this function takes the corresponding components of "ms.dVdq" and writes
// them into "m.dVdq" vectors for each monomer 
void extract_dVdq_and_write_into_monomers(MoleculeSystem *ms);

// 20.12.2024 NOTE: 
// the MoleculeSystem struct needs to be passed as void* into "rhs" function.
// Then, the first step is to write the phase point (N_Vector y) into the fields of MoleculeSystem
// using the method "put_qp_into_ms".
int rhs(realtype t, N_Vector y, N_Vector ydot, void *data);

Array compute_numerical_derivatives(double (*f)(double *q), double *q, size_t len, size_t order);
Array compute_numerical_rhs(MoleculeSystem *ms, size_t order); 
gsl_matrix* compute_numerical_jac(void (*transform_angles)(double *qlab, double *qmol), double *qlab, size_t ninput_coordinates, size_t noutput_coordinates, size_t order);

double kinetic_energy(MoleculeSystem *ms);
double Hamiltonian(MoleculeSystem *ms);

double generate_normal(double sigma); 

void q_generator(MoleculeSystem *ms, CalcParams *params);
void p_generator(MoleculeSystem *ms, double Temperature);
bool reject(MoleculeSystem *ms, double Temperature, double pesmin);

void j_monomer(Monomer m, double j[3]);
double torque_monomer(Monomer m); 
double find_closest_integer(double j);
double find_closest_half_integer(double j); 

void invert_momenta(MoleculeSystem *ms);

void calculate_M0(MoleculeSystem *ms, CalcParams *params, double Temperature, double *m, double *q);
void compute_dHdp(MoleculeSystem *ms, gsl_matrix* dHdp); 
void calculate_M2(MoleculeSystem *ms, CalcParams *params, double Temperature, double *m, double *q);
double analytic_full_partition_function_by_V(MoleculeSystem *ms, double Temperature);

void recv_histogram_and_append(Arena *a, int source, gsl_histogram **h);

#ifdef USE_MPI
void mpi_calculate_M0(MoleculeSystem *ms, CalcParams *params, double Temperature, double *m, double *q);
void mpi_calculate_M2(MoleculeSystem *ms, CalcParams *params, double Temperature, double *m, double *q);

CFnc      calculate_correlation_and_save(MoleculeSystem *ms, CalcParams *params, double Temperature);
SFnc      calculate_spectral_function_using_prmu_representation_and_save(MoleculeSystem *ms, CalcParams *params, double Temperature);
CFncArray calculate_correlation_array_and_save(MoleculeSystem *ms, CalcParams *params, double base_temperature); 
#endif // USE_MPI
    
double* linspace(double start, double end, size_t n);
double* arena_linspace(Arena *a, double start, double end, size_t n); 
size_t* linspace_size_t(Arena *a, size_t start, size_t end, size_t n);

double integrate_composite_simpson(double *x, double *y, size_t len); 
double compute_Mn_from_sf_using_classical_detailed_balance(SFnc sf, size_t n);
double compute_Mn_from_sf_using_quantum_detailed_balance(SFnc sf, size_t n);

void free_cfnc(CFnc cf);
void free_cfnc_array(CFncArray ca);
void free_sfnc(SFnc cf);
void free_spectrum(Spectrum sp); 

// ----------------------------------------------------------
// String manipulation
// ----------------------------------------------------------
/* 
 * This struct represents a resizable buffer designed to build strings 
 * dynamically. It stores characters in a contiguous block of memory, 
 * allowing for manipulation of strings.
 */
typedef struct {
    char *items;
    size_t count;
    size_t capacity;
} String_Builder;

void sb_append(String_Builder *sb, const char *line, size_t n);
void sb_append_cstring(String_Builder *sb, const char *line);
void sb_append_format(String_Builder *sb, const char *format, ...);
void sb_append_seconds_as_datetime_string(String_Builder *sb, int seconds);
void sb_reset(String_Builder *sb);
void sb_free(String_Builder *sb);
// ----------------------------------------------------------


// ----------------------------------------------------------
// Processing the results   
// ----------------------------------------------------------
int save_correlation_function(FILE *fp, CFnc cf);
void save_spectral_function(FILE *fp, SFnc sf, CalcParams *params);

bool read_correlation_function(const char *filename, String_Builder *sb, CFnc *cf); 
bool read_spectral_function(const char *filename, String_Builder *sb, SFnc *sf); 
bool writetxt(const char *filename, double *x, double *y, size_t len, const char *header); 
#define average_correlation_functions(average, ...) average_correlation_functions__impl(average, sizeof((CFnc[]){__VA_ARGS__}) / sizeof(CFnc), __VA_ARGS__)
int average_correlation_functions__impl(CFnc *average, int arg_count,  ...);

int assert_float_is_equal_to(double estimate, double true_value, double abs_tolerance);

extern WingParams INIT_WP;
double wingmodel(WingParams *wp, double t);

int wingmodel_f(const gsl_vector* x, void* data, gsl_vector* f);
int wingmodel_df(const gsl_vector* x, void* data, gsl_matrix * J);

void gsl_nonlinear_opt(size_t n, double* x, double* y, WingParams *wing_params);
WingParams fit_baseline(CFnc *cf, size_t EXT_RANGE_MIN);

void connes_apodization(Array a, double sampling_time); 

gsl_histogram* gsl_histogram_extend_right(gsl_histogram* h, int add_bins);

/* Uses bit hack taken from: http://www.graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2 */ 
static inline bool is_power_of_two(size_t n) {
    return n && !(n & (n - 1));
}

static inline unsigned int round_to_next_power_of_two(unsigned int n)
/*
 * rounds up to the next highest power of 2 using a clever bit hack
 * Taken from: 
 *  https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
 * 
 * works on 32-bit numbers
 */ 
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;

    return n;
}

double* pad_to_power_of_two(double* v, size_t len, size_t* padded_len); 

double *dct(double *v, size_t len);
double *idct(double *v, size_t len);

SFnc dct_numeric_sf(CFnc cf, WingParams *wp);
CFnc dct_sf_to_cf(SFnc sf);
SFnc idct_cf_to_sf(CFnc cf);

CFnc egelstaff_time_transform(CFnc cf, bool frommhold_renormalization); 

SFnc desymmetrize_d2(SFnc sf); 
SFnc desymmetrize_schofield(SFnc sf); 
SFnc desymmetrize_egelstaff(SFnc sf);
SFnc desymmetrize_frommhold(SFnc sf);

Spectrum compute_alpha(SFnc sf); 
    
#ifdef __cplusplus
}
#endif

#endif // HAWAII_H_
