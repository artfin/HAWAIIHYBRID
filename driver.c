#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <inttypes.h>
#include <dlfcn.h>

#define USE_MPI
#include "hawaii.h"
#include "hep_hawaii.h"

// TODO: investigate the UNREACHABLE
// &INPUT
//  PAIR_STATE = ALL
//  CALCULATION_TYPE = CORRELATION_SINGLE
//  TEMPERATURE = 300
//&END



// TODO: add 'where_begin' field for token start in the Lexer

// TODO: add type of calculation for running trajectory from specific phase-point 

// TODO: we may want to have some predefined constants like 'l_H2' or 'm_H2'
// and have a library of them (we actually have one in constants.h) instead of specifying them as floats 
// in the input file

typedef struct {
    double reduced_mass;
    double Temperature;
    double *temperatures;
    size_t num_temperatures;
    const char* so_potential;
    const char* so_dipole_1;
    const char* so_dipole_2;
} InputBlock;

typedef struct {
    const char *input_path;
    int line_number;
    int line_offset;
} Loc;


typedef enum {
    TOKEN_EOF,
    TOKEN_BLOCK,
    TOKEN_KEYWORD,
    TOKEN_STRING,
    TOKEN_DQSTRING,
    TOKEN_INTEGER,
    TOKEN_FLOAT,
    TOKEN_BOOLEAN,
    TOKEN_OCURLY,
    TOKEN_CCURLY,
    TOKEN_OPAREN,
    TOKEN_CPAREN,
    TOKEN_COMMA,
    TOKEN_EQ,
    TOKEN_FUNCALL,
    TOKEN_COUNT, // sentinel value for array size
} Token_Type;

const char *TOKEN_TYPES[TOKEN_COUNT] = {
      [TOKEN_EOF]          = "EOF",
      [TOKEN_BLOCK]        = "BLOCK",
      [TOKEN_KEYWORD]      = "KEYWORD",
      [TOKEN_STRING]       = "STRING",
      [TOKEN_DQSTRING]     = "DOUBLE-QUOTED STRING",
      [TOKEN_INTEGER]      = "INTEGER",
      [TOKEN_FLOAT]        = "FLOAT",
      [TOKEN_BOOLEAN]      = "BOOLEAN",
      [TOKEN_OCURLY]       = "{",
      [TOKEN_CCURLY]       = "}",
      [TOKEN_OPAREN]       = "(",
      [TOKEN_CPAREN]       = ")",
      [TOKEN_COMMA]        = ",",
      [TOKEN_EQ]           = "=",
      [TOKEN_FUNCALL]      = "FUNCTION CALL",
}; 
static_assert(TOKEN_COUNT == 15, "");

const char *PUNCTS[TOKEN_COUNT] = {
    [TOKEN_COMMA]  = ",",
    [TOKEN_OCURLY] = "{",
    [TOKEN_CCURLY] = "}",
    [TOKEN_OPAREN] = "(",
    [TOKEN_CPAREN] = ")",
    [TOKEN_EQ]     = "=",
};

const char *BLOCK_NAMES[] = {
    [0] = "&INPUT",
    [1] = "&MONOMER",
    [2] = "&PROCESSING",
    [3] = "&END",
};
static_assert(sizeof(BLOCK_NAMES)/sizeof(BLOCK_NAMES[0]) == 4, "");

const char *BOOLEAN_AS_STR[] = {
    [0] = "FALSE",
    [1] = "TRUE",
};
static_assert(sizeof(BOOLEAN_AS_STR)/sizeof(BOOLEAN_AS_STR[0]) == 2, "");

typedef enum {
    /* INPUT BLOCK */
    KEYWORD_CALCULATION_TYPE,
    KEYWORD_PAIR_STATE,
    KEYWORD_PAIR_REDUCED_MASS,
    KEYWORD_SO_POTENTIAL,
    KEYWORD_SO_DIPOLE,
    KEYWORD_SO_DIPOLE_1,
    KEYWORD_SO_DIPOLE_2,

    KEYWORD_TEMPERATURE,
    KEYWORD_TEMPERATURES,
    KEYWORD_SATELLITE_TEMPERATURES,
    
    KEYWORD_NITERATIONS,
    KEYWORD_TOTAL_TRAJECTORIES,
    
    KEYWORD_CVODE_TOLERANCE,
    KEYWORD_SAMPLING_TIME,
    KEYWORD_MAXTRAJECTORYLENGTH,
    KEYWORD_SAMPLER_RMIN,
    KEYWORD_SAMPLER_RMAX,
    KEYWORD_PESMIN,

    KEYWORD_INITIALM0_NPOINTS,
    KEYWORD_INITIALM2_NPOINTS,
    KEYWORD_HEP_M0_NITERATIONS,
    KEYWORD_HEP_M0_NPOINTS,
    KEYWORD_HEP_M2_NITERATIONS,
    KEYWORD_HEP_M2_NPOINTS,
    KEYWORD_HEP_PPF_NITERATIONS,
    KEYWORD_HEP_PPF_NPOINTS,
    
    KEYWORD_SF_FILENAME,
    KEYWORD_CF_FILENAME,
    KEYWORD_CF_FILENAMES,
    KEYWORD_R0,
    KEYWORD_RCUT,
    KEYWORD_PARTIAL_PARTITION_FUNCTION_RATIO,
    KEYWORD_PARTIAL_PARTITION_FUNCTION_RATIOS,
    KEYWORD_APPROXIMATEFREQUENCYMAX,
    KEYWORD_ODD_J_SPIN_WEIGHT,
    KEYWORD_EVEN_J_SPIN_WEIGHT,
    KEYWORD_USE_ZIMMERMANN_TRICK,
    KEYWORD_AVERAGE_TIME_BETWEEN_COLLISIONS,

    /* MONOMER BLOCK */
    KEYWORD_MONOMER_TYPE,
    KEYWORD_DJ,
    KEYWORD_II,

    KEYWORD_INITIAL_J,
    
    KEYWORD_TORQUE_CACHE_LEN,
    KEYWORD_TORQUE_LIMIT,
    
    KEYWORD_NSWITCH_HISTOGRAM_BINS, 
    KEYWORD_NSWITCH_HISTOGRAM_MAX,
    KEYWORD_NSWITCH_HISTOGRAM_FILENAME,
    
    KEYWORD_JINI_HISTOGRAM_BINS,
    KEYWORD_JINI_HISTOGRAM_MAX,
    KEYWORD_JINI_HISTOGRAM_FILENAME,
    
    KEYWORD_JFIN_HISTOGRAM_BINS,
    KEYWORD_JFIN_HISTOGRAM_MAX,
    KEYWORD_JFIN_HISTOGRAM_FILENAME,

    /* PROCESSING BLOCK */
    KEYWORD_SPECTRUM_FREQUENCY_MAX,

    KEYWORD_COUNT,
} Keyword;

const char* KEYWORDS[KEYWORD_COUNT] = {
    [KEYWORD_CALCULATION_TYPE]                = "CALCULATION_TYPE",
    [KEYWORD_PAIR_STATE]                      = "PAIR_STATE",
    [KEYWORD_PAIR_REDUCED_MASS]               = "PAIR_REDUCED_MASS",
    [KEYWORD_SO_POTENTIAL]                    = "SO_POTENTIAL",
    [KEYWORD_SO_DIPOLE]                       = "SO_DIPOLE",
    [KEYWORD_SO_DIPOLE_1]                     = "SO_DIPOLE.1",
    [KEYWORD_SO_DIPOLE_2]                     = "SO_DIPOLE.2",
    [KEYWORD_TEMPERATURE]                     = "TEMPERATURE",
    [KEYWORD_TEMPERATURES]                    = "TEMPERATURES",
    [KEYWORD_SATELLITE_TEMPERATURES]          = "SATELLITE_TEMPERATURES",
    [KEYWORD_NITERATIONS]                     = "NITERATIONS",
    [KEYWORD_TOTAL_TRAJECTORIES]              = "TOTAL_TRAJECTORIES",
    [KEYWORD_CVODE_TOLERANCE]                 = "CVODE_TOLERANCE",
    [KEYWORD_SAMPLING_TIME]                   = "SAMPLING_TIME",
    [KEYWORD_MAXTRAJECTORYLENGTH]             = "MAXTRAJECTORYLENGTH",
    [KEYWORD_SAMPLER_RMIN]                    = "SAMPLER_RMIN",
    [KEYWORD_SAMPLER_RMAX]                    = "SAMPLER_RMAX",
    [KEYWORD_PESMIN]                          = "PESMIN",
    [KEYWORD_INITIALM0_NPOINTS]               = "INITIALM0_NPOINTS",
    [KEYWORD_INITIALM2_NPOINTS]               = "INITIALM2_NPOINTS",
    [KEYWORD_HEP_M0_NPOINTS]                  = "HEP_M0_NPOINTS",
    [KEYWORD_HEP_M0_NITERATIONS]              = "HEP_M0_NITERATIONS",
    [KEYWORD_HEP_M2_NPOINTS]                  = "HEP_M2_NPOINTS",
    [KEYWORD_HEP_M2_NITERATIONS]              = "HEP_M2_NITERATIONS",
    [KEYWORD_HEP_PPF_NITERATIONS]             = "HEP_PPF_NITERATIONS",
    [KEYWORD_HEP_PPF_NPOINTS]                 = "HEP_PPF_NPOINTS",
    [KEYWORD_SF_FILENAME]                     = "SF_FILENAME",
    [KEYWORD_CF_FILENAME]                     = "CF_FILENAME",
    [KEYWORD_CF_FILENAMES]                    = "CF_FILENAMES",
    [KEYWORD_R0]                              = "R0",
    [KEYWORD_RCUT]                            = "RCUT",
    [KEYWORD_PARTIAL_PARTITION_FUNCTION_RATIO] = "PARTIAL_PARTITION_FUNCTION_RATIO",
    [KEYWORD_PARTIAL_PARTITION_FUNCTION_RATIOS] = "PARTIAL_PARTITION_FUNCTION_RATIOS",
    [KEYWORD_APPROXIMATEFREQUENCYMAX]         = "APPROXIMATEFREQUENCYMAX",
    [KEYWORD_ODD_J_SPIN_WEIGHT]               = "ODD_J_SPIN_WEIGHT",
    [KEYWORD_EVEN_J_SPIN_WEIGHT]              = "EVEN_J_SPIN_WEIGHT",
    [KEYWORD_USE_ZIMMERMANN_TRICK]            = "USE_ZIMMERMANN_TRICK",
    [KEYWORD_AVERAGE_TIME_BETWEEN_COLLISIONS] = "AVERAGE_TIME_BETWEEN_COLLISIONS",
    /* MONOMER BLOCK */
    [KEYWORD_MONOMER_TYPE]                    = "MONOMER_TYPE",
    [KEYWORD_DJ]                              = "DJ",
    [KEYWORD_II]                              = "II",
    [KEYWORD_INITIAL_J]                       = "INITIAL_J",
    [KEYWORD_TORQUE_CACHE_LEN]                = "TORQUE_CACHE_LEN",
    [KEYWORD_TORQUE_LIMIT]                    = "TORQUE_LIMIT",
    [KEYWORD_NSWITCH_HISTOGRAM_BINS]          = "NSWITCH_HISTOGRAM_BINS",
    [KEYWORD_NSWITCH_HISTOGRAM_MAX]           = "NSWITCH_HISTOGRAM_MAX",
    [KEYWORD_NSWITCH_HISTOGRAM_FILENAME]      = "NSWITCH_HISTOGRAM_FILENAME",
    [KEYWORD_JINI_HISTOGRAM_BINS]             = "JINI_HISTOGRAM_BINS",
    [KEYWORD_JINI_HISTOGRAM_MAX]              = "JINI_HISTOGRAM_MAX",
    [KEYWORD_JINI_HISTOGRAM_FILENAME]         = "JINI_HISTOGRAM_FILENAME",
    [KEYWORD_JFIN_HISTOGRAM_BINS]             = "JFIN_HISTOGRAM_BINS",
    [KEYWORD_JFIN_HISTOGRAM_MAX]              = "JFIN_HISTOGRAM_MAX",
    [KEYWORD_JFIN_HISTOGRAM_FILENAME]         = "JFIN_HISTOGRAM_FILENAME",
    /* PROCESSING BLOCK */
    [KEYWORD_SPECTRUM_FREQUENCY_MAX]          = "SPECTRUM_FREQUENCY_MAX",
}; 
static_assert(KEYWORD_COUNT == 54, "");

Token_Type EXPECT_TOKEN_AFTER_KEYWORD[KEYWORD_COUNT] = {
    [KEYWORD_CALCULATION_TYPE]                = TOKEN_STRING,
    [KEYWORD_PAIR_STATE]                      = TOKEN_STRING,
    [KEYWORD_PAIR_REDUCED_MASS]               = TOKEN_FLOAT,
    [KEYWORD_SO_POTENTIAL]                    = TOKEN_DQSTRING,
    [KEYWORD_SO_DIPOLE]                       = TOKEN_DQSTRING,
    [KEYWORD_SO_DIPOLE_1]                     = TOKEN_DQSTRING,
    [KEYWORD_SO_DIPOLE_2]                     = TOKEN_DQSTRING,
    [KEYWORD_TEMPERATURE]                     = TOKEN_FLOAT,
    [KEYWORD_TEMPERATURES]                    = TOKEN_OCURLY,
    [KEYWORD_SATELLITE_TEMPERATURES]          = TOKEN_OCURLY,
    [KEYWORD_NITERATIONS]                     = TOKEN_INTEGER,
    [KEYWORD_TOTAL_TRAJECTORIES]              = TOKEN_INTEGER,
    [KEYWORD_CVODE_TOLERANCE]                 = TOKEN_FLOAT,
    [KEYWORD_SAMPLING_TIME]                   = TOKEN_FLOAT,
    [KEYWORD_MAXTRAJECTORYLENGTH]             = TOKEN_INTEGER,
    [KEYWORD_SAMPLER_RMIN]                    = TOKEN_FLOAT,
    [KEYWORD_SAMPLER_RMAX]                    = TOKEN_FLOAT,
    [KEYWORD_PESMIN]                          = TOKEN_FLOAT,
    [KEYWORD_INITIALM0_NPOINTS]               = TOKEN_INTEGER,
    [KEYWORD_INITIALM2_NPOINTS]               = TOKEN_INTEGER,
    [KEYWORD_HEP_M0_NPOINTS]                  = TOKEN_INTEGER,
    [KEYWORD_HEP_M0_NITERATIONS]              = TOKEN_INTEGER,
    [KEYWORD_HEP_M2_NPOINTS]                  = TOKEN_INTEGER,
    [KEYWORD_HEP_M2_NITERATIONS]              = TOKEN_INTEGER,
    [KEYWORD_HEP_PPF_NPOINTS]                 = TOKEN_INTEGER,
    [KEYWORD_HEP_PPF_NITERATIONS]             = TOKEN_INTEGER,
    [KEYWORD_SF_FILENAME]                     = TOKEN_DQSTRING,
    [KEYWORD_CF_FILENAME]                     = TOKEN_DQSTRING,
    [KEYWORD_CF_FILENAMES]                    = TOKEN_OCURLY,
    [KEYWORD_R0]                              = TOKEN_FLOAT,
    [KEYWORD_RCUT]                            = TOKEN_FLOAT,
    [KEYWORD_PARTIAL_PARTITION_FUNCTION_RATIO] = TOKEN_FLOAT,
    [KEYWORD_PARTIAL_PARTITION_FUNCTION_RATIOS] = TOKEN_OCURLY,
    [KEYWORD_APPROXIMATEFREQUENCYMAX]         = TOKEN_FLOAT,
    [KEYWORD_ODD_J_SPIN_WEIGHT]               = TOKEN_FLOAT,
    [KEYWORD_EVEN_J_SPIN_WEIGHT]              = TOKEN_FLOAT,
    [KEYWORD_USE_ZIMMERMANN_TRICK]            = TOKEN_BOOLEAN,
    [KEYWORD_AVERAGE_TIME_BETWEEN_COLLISIONS] = TOKEN_FLOAT,
    /* MONOMER BLOCK */
    [KEYWORD_MONOMER_TYPE]                    = TOKEN_STRING,
    [KEYWORD_DJ]                              = TOKEN_FLOAT,
    [KEYWORD_II]                              = TOKEN_OCURLY,
    [KEYWORD_INITIAL_J]                       = TOKEN_FLOAT,
    [KEYWORD_TORQUE_CACHE_LEN]                = TOKEN_INTEGER,
    [KEYWORD_TORQUE_LIMIT]                    = TOKEN_FLOAT,
    [KEYWORD_NSWITCH_HISTOGRAM_BINS]          = TOKEN_INTEGER,
    [KEYWORD_NSWITCH_HISTOGRAM_MAX]           = TOKEN_FLOAT,
    [KEYWORD_NSWITCH_HISTOGRAM_FILENAME]      = TOKEN_DQSTRING,
    [KEYWORD_JINI_HISTOGRAM_BINS]             = TOKEN_INTEGER,
    [KEYWORD_JINI_HISTOGRAM_MAX]              = TOKEN_FLOAT,
    [KEYWORD_JINI_HISTOGRAM_FILENAME]         = TOKEN_DQSTRING,
    [KEYWORD_JFIN_HISTOGRAM_BINS]             = TOKEN_INTEGER,
    [KEYWORD_JFIN_HISTOGRAM_MAX]              = TOKEN_FLOAT,
    [KEYWORD_JFIN_HISTOGRAM_FILENAME]         = TOKEN_DQSTRING,
    /* PROCESSING BLOCK */
    [KEYWORD_SPECTRUM_FREQUENCY_MAX]          = TOKEN_FLOAT, 
};

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} Satellite_Temperatures;

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} Temperatures;

typedef struct {
    char **items;
    size_t count;
    size_t capacity;
} CF_Filenames;

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} Partial_Partition_Function_Ratios;

typedef struct {
    char *input_stream;
    char *parse_point;
    char *eof;
    Token_Type token_type;
    char *token_start;
    size_t token_len;
    String_Builder string_storage;
    bool boolean_value;
    int64_t int_number;
    double double_number;
    Keyword keyword_type;
    Loc loc;
} Lexer;

typedef struct {
    const char *name;
    const char *arg;
    Loc loc;
} Funcall;

typedef struct {
    Funcall *items;
    size_t count;
    size_t capacity;
} Funcalls;

static const char *AVAILABLE_FUNCS[] = {
    "READ_CF",
    "WRITE_CF",
    "WRITE_SF",
    "WRITE_SPECTRUM",
    "FIT_BASELINE", // TODO: not descriptive
    "COMPUTE_M0_CLASSICAL_DETAILED_BALANCE",
    "COMPUTE_M2_CLASSICAL_DETAILED_BALANCE",
    "CF_TO_SF",
    "AVERAGE_CFS",
    "ALPHA",
    "D3",
    "DUP",
    "INT3",
};

typedef struct {
    double spectrum_frequency_max; 
    Funcalls fs;
} Processing_Params;

typedef union { 
    CFnc cf;
    SFnc sf;
    Spectrum sp;
} Stack_Item; 
 
typedef enum {
    STACK_ITEM_CF = 0,
    STACK_ITEM_SF,
    STACK_ITEM_SPECTRUM,
} Stack_Item_Type;

const char *STACK_ITEM_TYPES[] = {
    [STACK_ITEM_CF] = "CF", 
    [STACK_ITEM_SF] = "SF",
    [STACK_ITEM_SPECTRUM] = "SPECTRUM",
};

typedef struct {
    Stack_Item item;
    Stack_Item_Type typ;
    Loc loc;
} Tagged_Stack_Item;

typedef struct {
    Tagged_Stack_Item *items; 
    size_t count; 
    size_t capacity;

    // additional elements
    WingParams wing_params;
} Processing_Stack;

Lexer lexer_new(const char *input_path, const char *input_stream, char *eof)
{
    Lexer l = {0};
    
    l.input_stream = (char*) input_stream;
    l.eof = eof;
    l.parse_point = (char*) input_stream;
    l.loc.line_number = 1;
    l.loc.line_offset = 1;
    l.loc.input_path = input_path;
    memset(&l.string_storage, 0, sizeof(l.string_storage));

    return l;
}

bool is_eof(Lexer *l) {
    return l->parse_point >= l->eof;
}

char peek_char(Lexer *l) {
    if (is_eof(l)) return '\0';
    return *l->parse_point;
}

void skip_char(Lexer *l)
{
    assert(!is_eof(l));
    
    char c = *l->parse_point;
    l->parse_point++;
    l->token_len++;

    l->loc.line_offset++;

    if (c == '\n') {
        l->loc.line_offset = 0;
        l->loc.line_number++;
    } 
}

void skip_whitespaces(Lexer *l) {
    char c;
    while ((c = peek_char(l)) != '\0') {
        if (isspace(c)) {
            skip_char(l);
        } else {
            break;
        }
    }
}

bool skip_prefix(Lexer *l, const char *prefix) {
    if (!prefix) return false;

    char *saved_point = (*l).parse_point;

    while (*prefix) {
        char c = peek_char(l);
        if (c == '\0') {
            (*l).parse_point = saved_point;
            return false;
        }
        if (c != *prefix) {
            (*l).parse_point = saved_point;
            return false;
        }
        skip_char(l);
        prefix++;
    }

    return true;
}

void skip_until(Lexer *l, const char *prefix) {
    while (!is_eof(l) && !skip_prefix(l, prefix)) {
        skip_char(l);
    }
}

bool is_identifier_begin(char c) {
    return isalpha(c) || c == '_';
}

bool is_identifier(char c) {
    return isalnum(c) || c == '_' || c == '.';
}

bool get_token(Lexer *l) {
    while (true) {
        skip_whitespaces(l);

        if (skip_prefix(l, "!")) {
            skip_until(l, "\n");
            continue; 
        }

        break;
    }

    char c = peek_char(l);
    l->token_start = l->parse_point;
    l->token_len = 0;

    if (c == '\0') {
        l->token_type = TOKEN_EOF;
        return false; 
    }

    for (size_t i = 0; i < sizeof(PUNCTS)/sizeof(PUNCTS[0]); ++i) {
        const char *prefix = PUNCTS[i];
        if (skip_prefix(l, prefix)) {
            l->token_type = (Token_Type) i;
            l->token_len = 1;
            break;
        }
    } 

    if (c == '&') { // block begin
        (*l).token_type = TOKEN_BLOCK;
        (*l).string_storage.count = 0;
        
        da_append(&l->string_storage, c);
        skip_char(l); // consume '&'

        for ( c = peek_char(l); c != '\0'; c = peek_char(l)) {
            if (is_identifier(c)) {
                da_append(&l->string_storage, c);
                skip_char(l);
            } else {
                break;
            }
        }

        da_append(&(*l).string_storage, 0);
    }

    if (is_identifier_begin(c)) {
        (*l).token_type = TOKEN_STRING;
        (*l).string_storage.count = 0;
     
        c = peek_char(l); 
        while (is_identifier(c)) { 
            da_append(&l->string_storage, c);
            skip_char(l);
            c = peek_char(l);
        }

        da_append(&l->string_storage, 0);
   
        // looking for a "(" to gather the information that we have a function call  
        {
            char *saved_parse_point = l->parse_point;
            Loc saved_loc = l->loc;
            size_t saved_token_len = l->token_len;     

            // here we may consume new line and completely mess up the location of what is 
            // actually just a string and not a function call 
            skip_whitespaces(l);
            c = peek_char(l);
            if (c == '(') {
                l->token_type = TOKEN_FUNCALL;
            } else {
                l->parse_point = saved_parse_point;
                l->loc         = saved_loc;
                l->token_len   = saved_token_len;
            }
        }

        if (strcasecmp(l->string_storage.items, "true") == 0) {
            l->token_type = TOKEN_BOOLEAN;
            l->boolean_value = true;
        } else if (strcasecmp(l->string_storage.items, "false") == 0) {
            l->token_type = TOKEN_BOOLEAN;
            l->boolean_value = false;
        }

        for (size_t i = 0; i < sizeof(KEYWORDS)/sizeof(KEYWORDS[0]); ++i) {
            if (strcasecmp(l->string_storage.items, KEYWORDS[i]) == 0) {
                l->token_type = TOKEN_KEYWORD;
                l->keyword_type = i; 
            }
        }
    }

    if ((isdigit(c) != 0) || c == '-' || c == '+') {
        bool is_float = false;
        l->int_number = 0;
        l->double_number = 0;

        int sign = 1;
        if (c == '+') {
            skip_char(l);
        } else if (c == '-') {
            sign = -1;
            skip_char(l); 
        }

        for (c = peek_char(l); c != '\0'; c = peek_char(l)) {
            // TODO: check for overflow
            if (isdigit(c)) {
                l->int_number = l->int_number * 10 + (int)(c - '0');
                skip_char(l);
            } else if (c == '.' || c == 'e' || c == 'E') {
                is_float = true;
                break;
            } else if (c == '_') {
                skip_char(l);
            } else {
                break;
            }
        }

        l->int_number = sign * l->int_number;

        if (!is_float) {        
            l->token_type = TOKEN_INTEGER;
        } else {
            l->token_type = TOKEN_FLOAT;
            l->double_number = (double) l->int_number;
            if (peek_char(l) == '.') {
                skip_char(l); // consume '.'
                double fraction = 0.0;
                double divisor = 1.0;

                for (c = peek_char(l); isdigit(c); c = peek_char(l)) {
                    fraction = fraction*10.0 + (int)(c - '0');
                    divisor *= 10.0;
                    skip_char(l);
                } 
            
                l->double_number += sign * fraction/divisor;
            } 
            
            if (peek_char(l) == 'e' || peek_char(l) == 'E') {
                skip_char(l); // consume 'e' or 'E'

                int exp_sign = 1;
                c = peek_char(l);

                if (c == '+') {
                    skip_char(l);
                } else if (c == '-') {
                    exp_sign = -1;
                    skip_char(l);
                }

                int exponent_value = 0;
                for (c = peek_char(l); isdigit(c); c = peek_char(l)) {
                    exponent_value = exponent_value*10 + (c - '0');
                    skip_char(l);
                }
                
                l->double_number *= pow(10, exp_sign*exponent_value);
            }
        }
    }

    if (c == '"') {
        l->token_type = TOKEN_DQSTRING;
        skip_char(l); // consume beginning '"'
        
        (*l).string_storage.count = 0;
       
        // NOTE: this allows only for one-line strings
        // do not see the use for multi-line strings for now 
        for (c = peek_char(l); (c != '\0') && (c != '"') && !isspace(c); c = peek_char(l)) {
            da_append(&(*l).string_storage, c);
            skip_char(l);
        }
        
        if (c != '"') {
            PRINT0("ERROR: %s:%d:%d: unfinished string literal\n", l->loc.input_path, l->loc.line_number, l->loc.line_offset);
            exit(1);
        }
        skip_char(l); // consume finishing '"'
        da_append(&(*l).string_storage, 0);
    }

    // c = peek_char(l);
    // printf("UNKNOWN symbol: %c\n", c); 
    // UNREACHABLE("get_token");

    return true; 
}


void print_input_block(InputBlock *input_block) {
    printf("Input Block:\n");
    printf("  reduced_mass = %.5e\n", input_block->reduced_mass);
    printf("  Temperature = %.2e\n", input_block->Temperature);
    
    printf("  num_temperatures = %zu\n", input_block->num_temperatures);
    printf("  temperatures = {");
    for (size_t i = 0; i < input_block->num_temperatures; ++i) {
        printf("%.2e", input_block->temperatures[i]);
        if (i < input_block->num_temperatures - 1) printf(", ");
    }
    printf("}\n");
    
    printf("  so_potential = %s\n", input_block->so_potential);
    printf("  so_dipole_1  = %s\n", input_block->so_dipole_1);
    printf("  so_dipole_2  = %s\n", input_block->so_dipole_2);
}

void print_monomer(Monomer *monomer) {
    printf("Monomer:\n");
    printf("  t = %s\n", display_monomer_type(monomer->t));
    printf("  I = {%.5e, %.5e, %.5e}\n", monomer->II[0], monomer->II[1], monomer->II[2]);
    printf("  DJ = %.5e\n", monomer->DJ);  
    printf("  torque_cache_len        = %zu\n", monomer->torque_cache_len);
    printf("  torque_limit            = %.5e\n", monomer->torque_limit);
    printf("  jini_histogram_bins     = %zu\n", monomer->jini_histogram_bins);
    printf("  jini_histogram_max      = %.5e\n", monomer->jini_histogram_max);
    printf("  jfin_histogram_bins     = %zu\n", monomer->jfin_histogram_bins);
    printf("  jfin_histogram_max      = %.5e\n", monomer->jfin_histogram_max);
}

void print_params(CalcParams *params) {
    printf("Params:\n");
    printf("  pair state = %s\n", PAIR_STATES[params->ps]);
    printf("  calculation type = %s\n", CALCULATION_TYPES[params->calculation_type]);
    printf("  --- sampling ---\n");
    printf("  sampler_Rmin = %.5e\n", params->sampler_Rmin);
    printf("  sampler_Rmax = %.5e\n", params->sampler_Rmax);
    printf("  pesmin       = %.5e\n", params->pesmin);
    printf("  --- initial spectral moments check ---\n");
    printf("  initialM0_npoints = %zu\n", params->initialM0_npoints);
    printf("  initialM2_npoints = %zu\n", params->initialM2_npoints);
    printf("  partial_partition_function_ratio = %.5e\n", params->partial_partition_function_ratio);
    printf("  --- weights to factor in spin statistics ---\n"); 
    printf("  odd_j_spin_weight  = %.5e\n", params->odd_j_spin_weight);
    printf("  even_j_spin_weight = %.5e\n", params->even_j_spin_weight);
    printf("  --- trajectory ---\n"); 
    printf("  sampling_time = %.5e\n", params->sampling_time);
    printf("  MaxTrajectoryLength = %zu\n", params->MaxTrajectoryLength);
    printf("  allow_truncating_trajectories_at_length_limit = %d\n", params->allow_truncating_trajectories_at_length_limit);
    printf("  cvode_tolerance = %.5e\n", params->cvode_tolerance);
    printf("  --- iteration parameters  ---\n"); 
    printf("  niterations = %zu\n", params->niterations);
    printf("  total_trajectories = %zu\n", params->total_trajectories); 
    printf("  --- correlation function and correlation function array calculation ONLY --- \n"); 
    printf("  cf_filename = %s\n", params->cf_filename);
    printf("  Rcut = %.5e\n", params->Rcut);
    printf("  use_zimmermann_trick = %d\n", params->use_zimmermann_trick);
    printf("  --- pr/mu spectral function calculation ONLY --- \n"); 
    printf("  sf_filename = %s\n", params->sf_filename);
    printf("  ApproximateFrequencyMax = %.5e\n", params->ApproximateFrequencyMax);
    printf("  R0 = %.5e\n", params->R0);
    printf("  average_time_between_collisions = %.5e\n", params->average_time_between_collisions);
    printf("  --- correlation function array calculation ONLY --- \n"); 

    printf("  num_satellite_temperatures = %zu\n", params->num_satellite_temperatures);
    printf("  satellite_temperatures = {");
    for (size_t i = 0; i < params->num_satellite_temperatures; ++i) {
        printf("%.2e", params->satellite_temperatures[i]);
        if (i < params->num_satellite_temperatures - 1) printf(", ");
    }
    printf("}\n");

    printf("  cf_filenames = {");
    for (size_t i = 0; i < params->num_satellite_temperatures; ++i) {
        printf("%s", params->cf_filenames[i]);
        if (i < params->num_satellite_temperatures - 1) printf(", ");
    }
    printf("}\n");
 
    printf("  partial_partition_function_ratios = {");
    for (size_t i = 0; i < params->num_satellite_temperatures; ++i) {
        printf("%.5e", params->partial_partition_function_ratios[i]);
        if (i < params->num_satellite_temperatures - 1) printf(", ");
    }    
    printf("}\n");
}

void print_processing_params(Processing_Params *processing_params) {
    printf("Processing_Params:\n");
    printf("  spectrum_frequency_max: %.5e\n", processing_params->spectrum_frequency_max);

    for (size_t i = 0; i < processing_params->fs.count; ++i) {
        printf("  funcall: %s(%s)\n", processing_params->fs.items[i].name, processing_params->fs.items[i].arg);
    }
}

bool read_entire_file(const char *path, String_Builder *sb)
{
    bool result = true;

    FILE *fp = fopen(path, "rb");
    if (fp == NULL)                 return_defer(false);
    if (fseek(fp, 0, SEEK_END) < 0) return_defer(false);
    long m = ftell(fp);
    
    if (m < 0)                      return_defer(false);
    if (fseek(fp, 0, SEEK_SET) < 0) return_defer(false);

    sb_reserve(sb, m + 1);
    fread(sb->items + sb->count, m, 1, fp);
    sb->count += m;
    sb->items[sb->count] = '\0';

    if (ferror(fp)) {
       return_defer(false);
    }

defer:
    if (!result) PRINT0("ERROR: Could not read file %s: %s", path, strerror(errno));
    if (fp) fclose(fp);
    return result;
}

void print_lexeme(Lexer *l) {
    switch (l->token_type) {
      case TOKEN_KEYWORD:  PRINT0("KEYWORD: %s\n", l->string_storage.items); break;
      case TOKEN_STRING:   PRINT0("%s: '%s'\n", TOKEN_TYPES[l->token_type], l->string_storage.items); break;
      case TOKEN_BLOCK:    PRINT0("%s: %s\n", TOKEN_TYPES[l->token_type], l->string_storage.items); break;
      case TOKEN_INTEGER:  PRINT0("%s: %"PRId64"\n", TOKEN_TYPES[l->token_type], l->int_number); break;
      case TOKEN_FLOAT:    PRINT0("%s: %.16e\n", TOKEN_TYPES[l->token_type], l->double_number); break;
      case TOKEN_DQSTRING: PRINT0("%s: \"%s\"\n", TOKEN_TYPES[l->token_type], l->string_storage.items); break;
      case TOKEN_BOOLEAN:  PRINT0("%s: %d\n", TOKEN_TYPES[l->token_type], l->boolean_value); break;
      case TOKEN_EOF:      PRINT0("EOF\n"); break;
      case TOKEN_OCURLY:   PRINT0("TOKEN: %s\n", TOKEN_TYPES[l->token_type]); break; 
      case TOKEN_CCURLY:   PRINT0("TOKEN: %s\n", TOKEN_TYPES[l->token_type]); break; 
      case TOKEN_COMMA:    PRINT0("TOKEN: %s\n", TOKEN_TYPES[l->token_type]); break; 
      case TOKEN_EQ:       PRINT0("TOKEN: %s\n", TOKEN_TYPES[l->token_type]); break; 
      case TOKEN_COUNT:    UNREACHABLE("print_lexeme");
      default: {
        PRINT0("TOKEN: %s\n", TOKEN_TYPES[l->token_type]);
        break;
      }
    }
}

void print_lexemes(Lexer *l)
{
    for (; get_token(l); ) {
       print_lexeme(l); 
    } 

    assert(get_token(l) == TOKEN_EOF);
}

void expect_one_of_tokens(Lexer *l, int count, ...) {
    Token_Type t = l->token_type;

    va_list args;
    va_start(args, count);
    
    String_Builder sb = {0};

    for (int i = 0; i < count; ++i) {
        Token_Type expected = va_arg(args, Token_Type);
        if (t == expected) {
            va_end(args);
            return;
        }
    }

    va_start(args, count); // reset va_args

    bool print_boolean_hint = false;

    for (int i = 0; i < count; ++i) {
        Token_Type expected = va_arg(args, Token_Type);
        if (expected == TOKEN_BOOLEAN) {
            print_boolean_hint = true;
        }

        sb_append_format(&sb, "'%s'", TOKEN_TYPES[expected]);

        if (i < count - 1) sb_append_cstring(&sb, " or "); 
    }
    
    va_end(args);

    PRINT0("ERROR: %s:%d:%d: expected one of [%s] but got '%s'\n", 
            l->loc.input_path, l->loc.line_number, l->loc.line_offset,
            sb.items, TOKEN_TYPES[t]);
        
    if (print_boolean_hint) {
        PRINT0("Use TRUE and FALSE for boolean values\n");
    }

    sb_free(&sb); 
    exit(1);
}

void expect_token(Lexer *l, Token_Type expected) 
{
    Token_Type t = l->token_type;

    if (t != expected) {
        switch (t) {
            case TOKEN_KEYWORD:
            case TOKEN_BLOCK:
            case TOKEN_STRING: {
               PRINT0("ERROR: %s:%d:%d: expected token of type %s but got '%s' of type %s\n", 
                        l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len + 1,
                        TOKEN_TYPES[expected], l->string_storage.items, TOKEN_TYPES[t]);
               break; 
           }
           case TOKEN_DQSTRING: {
               PRINT0("ERROR: %s:%d:%d: expected token of type %s but got \"%s\" of type %s\n", 
                        l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len + 1,
                        TOKEN_TYPES[expected], l->string_storage.items, TOKEN_TYPES[t]);
               break; 
           }
           case TOKEN_COMMA:
           case TOKEN_OCURLY:
           case TOKEN_CCURLY:
           case TOKEN_OPAREN:
           case TOKEN_CPAREN:
           case TOKEN_EQ: {
               PRINT0("ERROR: %s:%d:%d: expected token of type %s but got punctuation token \"%s\" \n", 
                        l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len + 1,
                        TOKEN_TYPES[expected], TOKEN_TYPES[t]);
               break; 
           }
           case TOKEN_EOF: {
               PRINT0("ERROR: %s:%d:%d: expected token of type %s but got %s\n", 
                        l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len + 1,
                        TOKEN_TYPES[expected], TOKEN_TYPES[t]);
               break; 
           }
           case TOKEN_BOOLEAN: {
               PRINT0("ERROR: %s:%d:%d: expected token of type %s but got %s of type %s\n", 
                        l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len + 1,
                        TOKEN_TYPES[expected], BOOLEAN_AS_STR[l->boolean_value], TOKEN_TYPES[t]);
               break; 
           }
           case TOKEN_FLOAT: {
               PRINT0("ERROR: %s:%d:%d: expected token of type %s but got %.5e of type %s\n", 
                        l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len + 1,
                        TOKEN_TYPES[expected], l->double_number, TOKEN_TYPES[t]);
               break; 
           }
           case TOKEN_FUNCALL: { 
                PRINT0("ERROR: %s:%d:%d: expected token of type %s but got %s of type %s\n",
                       l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len + 1,
                       TOKEN_TYPES[expected], l->string_storage.items, TOKEN_TYPES[t]);
                break;
           }
           case TOKEN_COUNT: UNREACHABLE("expect_token"); 
           default: { 
               UNREACHABLE("expect_token"); 
           }
        }

        if (expected == TOKEN_BOOLEAN) {
            PRINT0("Use TRUE and FALSE for boolean values\n");
        }

        exit(1);
    }
}

void get_and_expect_token(Lexer *l, Token_Type token) {
    get_token(l);
    expect_token(l, token);
}

void parse_input_block(Lexer *l, InputBlock *input_block, CalcParams *params) 
{
    Satellite_Temperatures st = {0};
    Temperatures temperatures = {0}; 
    CF_Filenames cf_filenames = {0};
    Partial_Partition_Function_Ratios ppfs = {0};

    while (true) {
        get_token(l);
        if ((l->token_type == TOKEN_BLOCK) && (strcasecmp(l->string_storage.items, "&END") == 0)) {
            //print_lexeme(l);
            break; 
        }

        expect_token(l, TOKEN_KEYWORD);
        Keyword keyword_type = l->keyword_type;

        get_and_expect_token(l, TOKEN_EQ);
        
        Token_Type expect_token = EXPECT_TOKEN_AFTER_KEYWORD[keyword_type]; 
        get_and_expect_token(l, expect_token);
        
        switch (keyword_type) {
            case KEYWORD_CALCULATION_TYPE: {
                if (strcasecmp(l->string_storage.items, "PR_MU") == 0) {
                    params->calculation_type = CALCULATION_PR_MU;
                } else if (strcasecmp(l->string_storage.items, "CORRELATION_SINGLE") == 0) {
                    params->calculation_type = CALCULATION_CORRELATION_SINGLE;
                } else if (strcasecmp(l->string_storage.items, "CORRELATION_ARRAY") == 0) {
                    params->calculation_type = CALCULATION_CORRELATION_ARRAY;
                } else if (strcasecmp(l->string_storage.items, "CALCULATE_PHASE_SPACE_M0") == 0) {
                    params->calculation_type = CALCULATION_PHASE_SPACE_M0;
                } else if (strcasecmp(l->string_storage.items, "CALCULATE_PHASE_SPACE_M2") == 0) {
                    params->calculation_type = CALCULATION_PHASE_SPACE_M2;
                } else if (strcasecmp(l->string_storage.items, "PROCESSING") == 0) {
                    params->calculation_type = CALCULATION_PROCESSING;
                } else {
                    PRINT0("ERROR: %s:%d:%d: unknown calculation type '%s'\n", 
                            l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len+1, l->string_storage.items);
                    PRINT0("Available calculation types:\n");
                    
                    for (size_t i = 1; i < sizeof(CALCULATION_TYPES)/sizeof(CALCULATION_TYPES[0]); ++i) {
                        PRINT0("  %s\n", CALCULATION_TYPES[i]);
                    }

                    exit(1);
                }

                break;
            }
            case KEYWORD_PAIR_STATE: {
                if (strcasecmp(l->string_storage.items, "FREE_AND_METASTABLE") == 0) {
                    params->ps = PAIR_STATE_FREE_AND_METASTABLE;
                } else if (strcasecmp(l->string_storage.items, "BOUND") == 0) {
                    params->ps = PAIR_STATE_BOUND;
                } else if (strcasecmp(l->string_storage.items, "ALL") == 0) {
                    params->ps = PAIR_STATE_ALL;
                } else {
                    PRINT0("ERROR: %s:%d:%d: unknown pair state '%s'\n", 
                            l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len+1, l->string_storage.items);
                    PRINT0("Available pair states:\n");

                    for (size_t i = 0; i < sizeof(PAIR_STATES)/sizeof(PAIR_STATES[0]); ++i) {
                        PRINT0("  %s\n", PAIR_STATES[i]);
                    }

                    exit(1);
                }

                break;
            }
            case KEYWORD_PAIR_REDUCED_MASS:       input_block->reduced_mass = l->double_number; break;
            case KEYWORD_SO_POTENTIAL:            input_block->so_potential = strdup(l->string_storage.items); break;
            case KEYWORD_SO_DIPOLE: {
                 input_block->so_dipole_1 = strdup(l->string_storage.items); 
                 input_block->so_dipole_2 = strdup(l->string_storage.items); 
                 break;
            }
            case KEYWORD_SO_DIPOLE_1: input_block->so_dipole_1 = strdup(l->string_storage.items); break; 
            case KEYWORD_SO_DIPOLE_2: input_block->so_dipole_2 = strdup(l->string_storage.items); break;

            case KEYWORD_TEMPERATURE: input_block->Temperature = l->double_number; break;
            case KEYWORD_TEMPERATURES: {
               while(true) {
                   get_and_expect_token(l, TOKEN_FLOAT);
                   da_append(&temperatures, l->double_number);
   
                   get_token(l); 
                   expect_one_of_tokens(l, 2, TOKEN_COMMA, TOKEN_CCURLY);
                   if (l->token_type == TOKEN_CCURLY) break;
               }

               input_block->temperatures = temperatures.items; 
               input_block->num_temperatures = temperatures.count;
               break;
            }
            case KEYWORD_SATELLITE_TEMPERATURES:  {
               while(true) {
                   get_and_expect_token(l, TOKEN_FLOAT);
                   da_append(&st, l->double_number);
   
                   get_token(l); 
                   expect_one_of_tokens(l, 2, TOKEN_COMMA, TOKEN_CCURLY);
                   if (l->token_type == TOKEN_CCURLY) break;
               }

               params->satellite_temperatures = st.items; 
               params->num_satellite_temperatures = st.count;
               break;
            }
            case KEYWORD_NITERATIONS:             params->niterations = l->int_number; break;
            case KEYWORD_TOTAL_TRAJECTORIES:      params->total_trajectories = l->int_number; break;
            case KEYWORD_CVODE_TOLERANCE:         params->cvode_tolerance = l->double_number; break;
            case KEYWORD_SAMPLING_TIME:           params->sampling_time = l->double_number; break;
            case KEYWORD_MAXTRAJECTORYLENGTH:     params->MaxTrajectoryLength = l->int_number; break;
            case KEYWORD_SAMPLER_RMIN:            params->sampler_Rmin = l->double_number; break;
            case KEYWORD_SAMPLER_RMAX:            params->sampler_Rmax = l->double_number; break;
            case KEYWORD_PESMIN:                  params->pesmin = l->double_number; break;
            case KEYWORD_INITIALM0_NPOINTS:       params->initialM0_npoints = l->int_number; break;
            case KEYWORD_INITIALM2_NPOINTS:       params->initialM2_npoints = l->int_number; break;
            case KEYWORD_HEP_M0_NITERATIONS:      params->hep_m0_niterations = l->int_number; break;
            case KEYWORD_HEP_M0_NPOINTS:          params->hep_m0_npoints = l->int_number; break;
            case KEYWORD_HEP_M2_NITERATIONS:      params->hep_m2_niterations = l->int_number; break;
            case KEYWORD_HEP_M2_NPOINTS:          params->hep_m2_npoints = l->int_number; break;
            case KEYWORD_HEP_PPF_NITERATIONS:     params->hep_ppf_niterations = l->int_number; break;
            case KEYWORD_HEP_PPF_NPOINTS:         params->hep_ppf_npoints = l->int_number; break;
            case KEYWORD_SF_FILENAME:             params->sf_filename = strdup(l->string_storage.items); break;
            case KEYWORD_CF_FILENAME:             params->cf_filename = strdup(l->string_storage.items); break;
            case KEYWORD_CF_FILENAMES: {
                while (true) {
                    get_and_expect_token(l, TOKEN_DQSTRING);
                    da_append(&cf_filenames, NULL);
                    cf_filenames.items[cf_filenames.count-1] = strdup(l->string_storage.items);

                    get_token(l);
                    expect_one_of_tokens(l, 2, TOKEN_COMMA, TOKEN_CCURLY);
                    if (l->token_type == TOKEN_CCURLY) break;
                } 

                params->cf_filenames = (const char**) cf_filenames.items;
                break;
            }
            case KEYWORD_R0:                      params->R0 = l->double_number; break;
            case KEYWORD_RCUT:                    params->Rcut = l->double_number; break;
            case KEYWORD_PARTIAL_PARTITION_FUNCTION_RATIO: params->partial_partition_function_ratio = l->double_number; break;
            case KEYWORD_PARTIAL_PARTITION_FUNCTION_RATIOS: {
                while (true) {
                   get_and_expect_token(l, TOKEN_FLOAT);
                   da_append(&ppfs, l->double_number);
   
                   get_token(l); 
                   expect_one_of_tokens(l, 2, TOKEN_COMMA, TOKEN_CCURLY);
                   if (l->token_type == TOKEN_CCURLY) break;
                } 

                params->partial_partition_function_ratios = ppfs.items;
                break;
            }
            case KEYWORD_APPROXIMATEFREQUENCYMAX: params->ApproximateFrequencyMax = l->double_number; break;
            case KEYWORD_ODD_J_SPIN_WEIGHT:       params->odd_j_spin_weight = l->double_number; break;
            case KEYWORD_EVEN_J_SPIN_WEIGHT:      params->even_j_spin_weight = l->double_number; break;
            case KEYWORD_USE_ZIMMERMANN_TRICK:    params->use_zimmermann_trick = l->boolean_value; break;
            case KEYWORD_AVERAGE_TIME_BETWEEN_COLLISIONS: params->average_time_between_collisions = l->double_number; break;
            default: {
              PRINT0("ERROR: %s:%d:%d: keyword '%s' cannot be used within &INPUT block\n",
                     l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len+1, l->string_storage.items);
              exit(1);
            } 
        }
    }

    if (st.count != cf_filenames.count) {
        PRINT0("ERROR: Satellite temperature count (%s = %zu) does not match correlation filename count (%s = %zu)\n",
                KEYWORDS[KEYWORD_SATELLITE_TEMPERATURES], st.count,
                KEYWORDS[KEYWORD_CF_FILENAMES], cf_filenames.count);
        exit(1);
    }

}
            
void parse_monomer_block(Lexer *l, Monomer *m) 
{
    m->initial_j = -1.0;

    while (true) {
        get_token(l);
        if ((l->token_type == TOKEN_BLOCK) && (strcasecmp(l->string_storage.items, "&END") == 0)) {
            //print_lexeme(l);
            return;
        }

        expect_token(l, TOKEN_KEYWORD);
        Keyword keyword_type = l->keyword_type;
        Token_Type expect_token = EXPECT_TOKEN_AFTER_KEYWORD[keyword_type]; 

        get_and_expect_token(l, TOKEN_EQ);
        get_and_expect_token(l, expect_token);

        switch (keyword_type) {
            case KEYWORD_MONOMER_TYPE: {
                if (strcasecmp(l->string_storage.items, "ATOM") == 0) {
                    m->t = ATOM;
                } else if (strcasecmp(l->string_storage.items, "LINEAR_MOLECULE") == 0) {
                    m->t = LINEAR_MOLECULE;
                } else if (strcasecmp(l->string_storage.items, "LINEAR_MOLECULE_REQ_INTEGER") == 0) {
                    m->t = LINEAR_MOLECULE_REQ_INTEGER;
                } else if (strcasecmp(l->string_storage.items, "LINEAR_MOLECULE_REQ_HALFINTEGER") == 0) {
                    m->t = LINEAR_MOLECULE_REQ_HALFINTEGER;
                } else if (strcasecmp(l->string_storage.items, "ROTOR") == 0) {
                    m->t = ROTOR;
                } else {
                    PRINT0("ERROR: %s:%d:%d: unknown monomer type '%s'\n", 
                            l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len+1, l->string_storage.items);
                    PRINT0("Available monomer_block types:\n");
                    for (size_t i = 0; i < sizeof(MONOMER_TYPES)/sizeof(MONOMER_TYPES[0]); ++i) {
                        PRINT0(" %s\n", display_monomer_type(MONOMER_TYPES[i])); 
                    }
                    exit(1);
                }
                
                break;
            }
            case KEYWORD_DJ: m->DJ = l->double_number; break; 
            case KEYWORD_II: {
               get_and_expect_token(l, TOKEN_FLOAT);
               m->II[0] = l->double_number;
               get_and_expect_token(l, TOKEN_COMMA);
               
               get_and_expect_token(l, TOKEN_FLOAT);
               m->II[1] = l->double_number;
               get_and_expect_token(l, TOKEN_COMMA);
               
               get_and_expect_token(l, TOKEN_FLOAT);
               m->II[2] = l->double_number;

               get_and_expect_token(l, TOKEN_CCURLY);
               break;
            } 
            case KEYWORD_INITIAL_J:                  m->initial_j = l->double_number; break;
            case KEYWORD_TORQUE_CACHE_LEN:           m->torque_cache_len = l->int_number; break;
            case KEYWORD_TORQUE_LIMIT:               m->torque_limit = l->double_number; break;
            case KEYWORD_NSWITCH_HISTOGRAM_BINS:     m->nswitch_histogram_bins = l->int_number; break; 
            case KEYWORD_NSWITCH_HISTOGRAM_MAX:      m->nswitch_histogram_max = l->double_number; break; 
            case KEYWORD_NSWITCH_HISTOGRAM_FILENAME: m->nswitch_histogram_filename = strdup(l->string_storage.items); break;
            case KEYWORD_JINI_HISTOGRAM_BINS:        m->jini_histogram_bins = l->int_number; break;
            case KEYWORD_JINI_HISTOGRAM_MAX:         m->jini_histogram_max = l->double_number; break;
            case KEYWORD_JINI_HISTOGRAM_FILENAME:    m->jini_histogram_filename = strdup(l->string_storage.items); break;
            case KEYWORD_JFIN_HISTOGRAM_BINS:        m->jfin_histogram_bins = l->int_number; break;
            case KEYWORD_JFIN_HISTOGRAM_MAX:         m->jfin_histogram_max = l->double_number; break;
            case KEYWORD_JFIN_HISTOGRAM_FILENAME:    m->jfin_histogram_filename = strdup(l->string_storage.items); break;
 
            default: {
              PRINT0("ERROR: %s:%d:%d: keyword '%s' cannot be used within &MONOMER block\n",
                     l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len+1, l->string_storage.items);
              exit(1);
            } 
        }
    }
} 


void parse_processing_block(Lexer *l, Processing_Params *processing_params) 
{
    while (true) {
        get_token(l);
        if ((l->token_type == TOKEN_BLOCK) && (strcasecmp(l->string_storage.items, "&END") == 0)) {
            //print_lexeme(l);
            return;
        }

        expect_one_of_tokens(l, 3, TOKEN_KEYWORD, TOKEN_FUNCALL);

        switch (l->token_type) {
            case TOKEN_KEYWORD: {
                get_and_expect_token(l, TOKEN_EQ);
                get_and_expect_token(l, EXPECT_TOKEN_AFTER_KEYWORD[l->keyword_type]);

                switch (l->keyword_type) {
                    case KEYWORD_SPECTRUM_FREQUENCY_MAX: { 
                        processing_params->spectrum_frequency_max = l->double_number; break; 
                    }

                    default: {
                        PRINT0("ERROR: %s:%d:%d: keyword '%s' cannot be used within &PROCESSING block\n",
                               l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len+1, l->string_storage.items);
                        exit(1);
                    }
                }

                break;
            }
            case TOKEN_FUNCALL: {
                Funcall f = {0};

                bool funcall_is_available = false;
                for (size_t i = 0; i < sizeof(AVAILABLE_FUNCS)/sizeof(AVAILABLE_FUNCS[0]); ++i) {
                    if (strcasecmp(l->string_storage.items, AVAILABLE_FUNCS[i]) == 0) {
                        funcall_is_available = true;
                        break;
                    }
                }
                if (!funcall_is_available) {
                    PRINT0("ERROR: %s:%d:%d: funcall '%s' is unknown and cannot be used within PROCESSING block\n",
                            l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len+1, l->string_storage.items);

                    exit(1); 
                }

                f.name = strdup(l->string_storage.items);
                f.loc = (Loc) {
                    .input_path = strdup(l->loc.input_path),
                    .line_number = l->loc.line_number,
                    .line_offset = l->loc.line_offset-(int)l->token_len+1,
                }; 

                get_and_expect_token(l, TOKEN_OPAREN);
                
                get_token(l);
                expect_one_of_tokens(l, 2, TOKEN_CPAREN, TOKEN_DQSTRING);

                if (l->token_type == TOKEN_DQSTRING) {
                    f.arg = strdup(l->string_storage.items);
                    get_and_expect_token(l, TOKEN_CPAREN);
                }

                da_append(&processing_params->fs, f);
                break;
            }

            default: UNREACHABLE("parse_processing_block");
        } 
    }
} 

void parse_params(Lexer *l, CalcParams *calc_params, InputBlock *input_block, Monomer *m1, Monomer *m2, Processing_Params *processing_params)
{
    size_t monomer_blocks_count = 0;
    Monomer *m = m1;

    while (true) { 
        get_token(l);
        if (is_eof(l)) break;

        expect_token(l, TOKEN_BLOCK);

        if (strcasecmp(l->string_storage.items, "&END") == 0) {
            PRINT0("ERROR: %s:%d:%d: found '%s' without corresponding block beginning\n", 
                    l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len+1, l->string_storage.items);
            exit(1);
        } else if (strcasecmp(l->string_storage.items, "&INPUT") == 0) {
            parse_input_block(l, input_block, calc_params);
        } else if (strcasecmp(l->string_storage.items, "&MONOMER") == 0) {
            if (monomer_blocks_count >= 2) {
                PRINT0("ERROR: %s:%d:%d: only two &MONOMER blocks could be defined\n", l->loc.input_path, l->loc.line_number, l->loc.line_offset);
                exit(1);
            }

            parse_monomer_block(l, m);

            m = m2;
            monomer_blocks_count++; 
        } else if (strcasecmp(l->string_storage.items, "&PROCESSING") == 0) {
            parse_processing_block(l, processing_params);
        } else {
            PRINT0("ERROR: %s:%d:%d: found unknown block name '%s'\n",
                   l->loc.input_path, l->loc.line_number, l->loc.line_offset-(int)l->token_len+1,
                   l->string_storage.items);
            PRINT0("  Expected block names:\n");
            for (size_t i = 0; i < sizeof(BLOCK_NAMES)/sizeof(BLOCK_NAMES[0]); ++i) {
                PRINT0("    %s\n", BLOCK_NAMES[i]);
            }

            exit(1);
        }
    }
    
    if (calc_params->calculation_type == CALCULATION_NONE) {
        PRINT0("ERROR: Required field missing: '%s'\n", KEYWORDS[KEYWORD_CALCULATION_TYPE]);
        exit(1); 
    } else if ((calc_params->calculation_type != CALCULATION_PROCESSING) && 
               (calc_params->calculation_type != CALCULATION_PHASE_SPACE_M0)) {
        if (calc_params->ps == PAIR_STATE_NONE) {
            PRINT0("ERROR: Required field missing: '%s'\n", KEYWORDS[KEYWORD_PAIR_STATE]);
            exit(1);
        }
        
        if (input_block->Temperature <= 0.0) {
            PRINT0("ERROR: Required field missing: '%s'\n", KEYWORDS[KEYWORD_TEMPERATURE]);
            exit(1);
        }
    }
}
   
void stack_push(Processing_Stack *stack, Tagged_Stack_Item tagged_item) {
    da_append(stack, tagged_item);
}

void stack_push_with_type(Processing_Stack *stack, void *item, Stack_Item_Type typ, Loc *loc) 
{
    Tagged_Stack_Item tagged_item = {
        .typ = typ,
        .loc = *loc, 
    };

    switch (typ) {
        case STACK_ITEM_CF: memcpy(&tagged_item.item.cf, item, sizeof(CFnc)); break;
        case STACK_ITEM_SF: memcpy(&tagged_item.item.sf, item, sizeof(SFnc)); break;
        case STACK_ITEM_SPECTRUM: memcpy(&tagged_item.item.sp, item, sizeof(Spectrum)); break; 
    }
    
    da_append(stack, tagged_item); 
}

Tagged_Stack_Item stack_pop_with_type(Processing_Stack *stack, Loc loc) {
    if (stack->count < 1) {
        PRINT0("ERROR: %s:%d:%d: Cannot pop from an empty stack\n", 
                loc.input_path, loc.line_number, loc.line_offset);
        exit(1);
    }

    return stack->items[--stack->count];
}

Tagged_Stack_Item *stack_peek_with_type(Processing_Stack *stack, size_t i, Loc loc) {
    if (i >= stack->count) {
        PRINT0("ERROR: %s:%d:%d: Cannot peek element %zu on a stack with %zu elements\n", 
                loc.input_path, loc.line_number, loc.line_offset, i, stack->count);
        exit(1);
    }

    return &stack->items[stack->count-1 - i];
}

Tagged_Stack_Item *stack_peek_top_with_type(Processing_Stack *stack) {
    return &stack->items[stack->count - 1];
}

void expect_item_on_stack(Loc *pc_loc, Tagged_Stack_Item *tagged_item, Stack_Item_Type expected_type) {
    if (tagged_item->typ != expected_type) {
        PRINT0("ERROR: %s:%d:%d: corrupted stack: expected to find %s but found %s, which is created at\n", 
               pc_loc->input_path, pc_loc->line_number, pc_loc->line_offset, 
               STACK_ITEM_TYPES[expected_type], STACK_ITEM_TYPES[tagged_item->typ]);
        PRINT0("       %s:%d:%d\n", 
                tagged_item->loc.input_path, tagged_item->loc.line_number, tagged_item->loc.line_offset);
        exit(1);
    }
}

void expect_one_of_items_on_stack(Loc *pc_loc, Tagged_Stack_Item *tagged_item, int count, ...)
{
    va_list args;
    va_start(args, count);

    for (int i = 0; i < count; ++i) {
        Stack_Item_Type expected_typ = va_arg(args, Stack_Item_Type);
        if (tagged_item->typ == expected_typ) {
            va_end(args);
            return;
        }
    }
    
    // constructing error message 
    va_start(args, count); // reset va_args

    String_Builder sb_expected_types = {0};

    for (int i = 0; i < count; ++i) {
        Stack_Item_Type expected_typ = va_arg(args, Stack_Item_Type);

        sb_append_format(&sb_expected_types, "'%s'", STACK_ITEM_TYPES[expected_typ]);
        if (i < count - 1) sb_append_cstring(&sb_expected_types, " or ");
    }

    va_end(args);

    PRINT0("ERROR: %s:%d:%d: corrupted stack: expected to find one of [%s] but found '%s', which is created at\n",
            pc_loc->input_path, pc_loc->line_number, pc_loc->line_offset,
            sb_expected_types.items, STACK_ITEM_TYPES[tagged_item->typ]);
    PRINT0("       %s:%d:%d\n", tagged_item->loc.input_path, tagged_item->loc.line_number, tagged_item->loc.line_offset);

    sb_free(&sb_expected_types);
    exit(1);
}

void expect_arg(Funcall *func) {
    if (func->arg == NULL) {
        PRINT0("ERROR: %s:%d:%d: missing argument for %s\n",
                func->loc.input_path, func->loc.line_number, func->loc.line_offset, func->name);
        exit(1); 
    }
}

bool run_processing(Processing_Params *processing_params) {
    print_processing_params(processing_params);
    PRINT0("\n\n");

    Processing_Stack stack = {0};

    for (size_t pc = 0; pc < processing_params->fs.count; ++pc) {
        Funcall *func = &processing_params->fs.items[pc];
        const char *funcname = func->name;
        Loc *pc_loc = &func->loc;

        PRINT0("\n");
        PRINT0("[%zu] %s\n", pc+1, funcname);
        _print0_margin = 2;

        if (strcasecmp(funcname, "READ_CF") == 0) {
            expect_arg(func);
            const char *filename = processing_params->fs.items[pc].arg;

            CFnc cf = {0}; 
            
            if (!read_correlation_function(filename, NULL, &cf)) {
                PRINT0("ERROR: could not read the file '%s'!\n", filename);
                return false; 
            }

            //for (size_t i = 0; i < cf.len; ++i) {
            //    cf.t[i] = cf.t[i] * ATU; 
            //}

            stack_push_with_type(&stack, (void*) &cf, STACK_ITEM_CF, pc_loc);
       
        } else if (strcasecmp(funcname, "DUP") == 0) {
            Tagged_Stack_Item *tagged_item = stack_peek_top_with_type(&stack);
            INFO("Dupicating %s on processing stack\n", STACK_ITEM_TYPES[tagged_item->typ]);

            switch (tagged_item->typ) {
                case STACK_ITEM_CF: {
                    CFnc cf_copy = copy_cfnc(tagged_item->item.cf);
                    stack_push(&stack, (Tagged_Stack_Item) {
                        .item.cf = cf_copy,
                        .typ = STACK_ITEM_CF,
                        .loc = *pc_loc,
                    });
                    break;
                }
                case STACK_ITEM_SF: {
                    SFnc sf_copy = copy_sfnc(tagged_item->item.sf);
                    stack_push(&stack, (Tagged_Stack_Item) {
                        .item.sf = sf_copy,
                        .typ = STACK_ITEM_SF,
                        .loc = *pc_loc,
                    });
                    break;
                } 
                case STACK_ITEM_SPECTRUM: {
                    Spectrum sp_copy = copy_spectrum(tagged_item->item.sp);
                    stack_push(&stack, (Tagged_Stack_Item) {
                        .item.sp = sp_copy,
                        .typ = STACK_ITEM_SPECTRUM,
                        .loc = *pc_loc,
                    });
                    break;
                } 
            }
        } else if (strcasecmp(funcname, "INT3") == 0) {
            PRINT0("BREAKPOINT INTERRUPT ISSUED\n");
            PRINT0("  Stack trace:\n");
            for (size_t i = 0; i < stack.count; ++i) {
                Loc *loc = &stack.items[i].loc;
                PRINT0("    %zu: %s created at %s:%d:%d\n", 
                        i, STACK_ITEM_TYPES[stack.items[i].typ], loc->input_path, loc->line_number, loc->line_offset);
            }

            exit(1);  

        } else if (strcasecmp(funcname, "AVERAGE_CFS") == 0) {
            CFncs cfncs = {0};

            // peeking the CFnc elements on the stack to have their valid address in memory
            // to store these addresses in dynamic memory  
            for (size_t i = 0; i < stack.count; ++i) {
                Tagged_Stack_Item *tagged_item = stack_peek_with_type(&stack, i, *pc_loc);
                expect_item_on_stack(pc_loc, tagged_item, STACK_ITEM_CF);

                da_append(&cfncs, &tagged_item->item.cf);
                
                INFO("Popped CF (i: %zu) from processing stack: ntraj = %.4e\n", cfncs.count, da_last(&cfncs)->ntraj);
            }

            for (size_t i = 0; i < cfncs.count; ++i) {
                PRINT0("CF(%zu): %p => ntraj = %.4e\n", i, cfncs.items[i], cfncs.items[i]->ntraj);
            } 

            CFnc average = {0};
            if (!average_correlation_functions(&average, cfncs)) {
                PRINT0("ERROR: %s:%d:%d: An error occured during averaging of correlation functions\n",
                        pc_loc->input_path, pc_loc->line_number, pc_loc->line_offset);
                exit(1);
            }

            for (size_t i = 0; i < cfncs.count; ++i) {
                free_cfnc(*cfncs.items[i]);
            }
            stack.count = 0; // manually resetting the state of the stack instead of 'pop'

            stack_push_with_type(&stack, (void*) &average, STACK_ITEM_CF, pc_loc);

        } else if (strcasecmp(funcname, "COMPUTE_M0_CLASSICAL_DETAILED_BALANCE") == 0) {
            Tagged_Stack_Item tagged_item = stack_pop_with_type(&stack, *pc_loc);
            expect_one_of_items_on_stack(pc_loc, &tagged_item, 2, STACK_ITEM_CF, STACK_ITEM_SF);

            if (tagged_item.typ == STACK_ITEM_SF) {
                SFnc *sf = &tagged_item.item.sf;

                if (processing_params->spectrum_frequency_max > 0) {
                    assert(sf->len > 1);
                    double dnu = sf->nu[1] - sf->nu[0];
                    size_t n = (size_t) (processing_params->spectrum_frequency_max / dnu);
                    sf->len = (sf->len > n) ? n : sf->len;
                    
                    INFO("Truncating spectral function at max frequency = %.4e cm-1 (resulting npoints: %zu)\n", n*dnu, n);
                }

                double M0 = compute_Mn_from_sf_using_classical_detailed_balance(*sf, 0);
                INFO("M0 BASED ON SPECTRAL FUNCTION = %.5e\n", M0); 
            } else if (tagged_item.typ == STACK_ITEM_CF) {
                CFnc *cf = &tagged_item.item.cf;

                double M0;
                compute_Mn_from_cf_using_classical_detailed_balance(*cf, 0, &M0);
                INFO("M0 BASED ON CORRELATION FUNCTION = %.5e\n", M0);
            }
        
        } else if (strcasecmp(funcname, "COMPUTE_M2_CLASSICAL_DETAILED_BALANCE") == 0) {
            Tagged_Stack_Item tagged_item = stack_pop_with_type(&stack, *pc_loc);
            expect_one_of_items_on_stack(pc_loc, &tagged_item, 2, STACK_ITEM_CF, STACK_ITEM_SF);
            
            if (tagged_item.typ == STACK_ITEM_SF) { 
                SFnc *sf = &tagged_item.item.sf;

                if (processing_params->spectrum_frequency_max > 0) {
                    assert(sf->len > 1);
                    double dnu = sf->nu[1] - sf->nu[0];
                    size_t n = (size_t) (processing_params->spectrum_frequency_max / dnu);
                    sf->len = (sf->len > n) ? n : sf->len;
                    
                    INFO("Truncating spectral function at max frequency = %.4e cm-1 (resulting npoints: %zu)\n", n*dnu, n);
                }

                double M2 = compute_Mn_from_sf_using_classical_detailed_balance(*sf, 2);
                INFO("M2 = %.5e\n", M2); 
            } else if (tagged_item.typ == STACK_ITEM_CF) {
                CFnc *cf = &tagged_item.item.cf;
      
                double M2; 
                if (compute_Mn_from_cf_using_classical_detailed_balance(*cf, 2, &M2)) {
                    INFO("M2 BASED ON CORRELATION FUNCTION = %.5e\n", M2);
                }
            }

        } else if (strcasecmp(funcname, "FIT_BASELINE") == 0) {
            Tagged_Stack_Item tagged_item = stack_pop_with_type(&stack, *pc_loc);
            expect_item_on_stack(pc_loc, &tagged_item, STACK_ITEM_CF); 

            // TODO: make this parameter adjustable from &PROCESSING section          
            //size_t EXT_RANGE_MIN = 8192;
            size_t EXT_RANGE_MIN = 16384; 
            INFO("Using points starting from %zu for CF extrapolation\n", EXT_RANGE_MIN); 
            CFnc *cf = &tagged_item.item.cf; 

            // I don't want to create an additional Stack_Item for these parameters, so
            // just saving them into Processing_Stack to make them available in the 
            // 'CF_TO_SF' call 
            stack.wing_params = fit_baseline(cf, EXT_RANGE_MIN);

            for (size_t i = 0; i < cf->len; ++i) {
                double wing_value = wingmodel(&stack.wing_params, cf->t[i]); 
                //printf("i = %zu => CF = %.5e, wing = %.5e\n", i, cf->data[i], wing_value);
                cf->data[i] -= wing_value; 
            }

            stack_push_with_type(&stack, (void*) cf, STACK_ITEM_CF, pc_loc);

        } else if (strcasecmp(funcname, "CF_TO_SF") == 0) {
            Tagged_Stack_Item tagged_item = stack_pop_with_type(&stack, *pc_loc);
            expect_item_on_stack(pc_loc, &tagged_item, STACK_ITEM_CF); 

            CFnc *cf = &tagged_item.item.cf;
            SFnc sf = idct_cf_to_sf(*cf);
            free_cfnc(*cf);

            if ((stack.wing_params.A != 0) || (stack.wing_params.B != 0) || (stack.wing_params.C != 0)) {
                INFO("Adding Fourier image of the baseline to the numerical result\n");
                
                double Xscale = 1.0/LightSpeed_cm / ATU / 2.0 / M_PI;
                double Yscale = ATU*ADIPMOMU*ADIPMOMU / (4.0*M_PI*EPSILON0);
                
                stack.wing_params.A = stack.wing_params.A * Yscale * Xscale; 
                stack.wing_params.B = stack.wing_params.B * Xscale; 
                stack.wing_params.C = stack.wing_params.C * Yscale * Xscale;
             
                INFO("Adding the delta function as an exponential Q-branch\n"); 
                for (size_t i = 0; i < sf.len; ++i) { 
                    double wingmodel_image_val = (stack.wing_params.C+stack.wing_params.A)/stack.wing_params.A*wingmodel_image(&stack.wing_params, sf.nu[i]);
                    sf.data[i] += wingmodel_image_val;
                    if (i < 10) printf("nu = %.5e => sf = %.5e, wmi = %.5e\n", sf.nu[i], sf.data[i], wingmodel_image_val);
                }
            }

            stack_push_with_type(&stack, (void*) &sf, STACK_ITEM_SF, pc_loc);

        } else if (strcasecmp(funcname, "ALPHA") == 0) {
            Tagged_Stack_Item tagged_item = stack_pop_with_type(&stack, *pc_loc);
            expect_item_on_stack(pc_loc, &tagged_item, STACK_ITEM_SF);

            SFnc *sf = &tagged_item.item.sf;
            Spectrum sp = compute_alpha(*sf);
            free_sfnc(*sf);

            stack_push_with_type(&stack, (void*) &sp, STACK_ITEM_SPECTRUM, pc_loc);

        } else if (strcasecmp(funcname, "D3") == 0) {
            // TODO: allow expecting SF or Spectrum
            // TODO: implement another 'desymmetrize_schofield' that yields SF
            Tagged_Stack_Item tagged_item = stack_pop_with_type(&stack, *pc_loc);
            expect_item_on_stack(pc_loc, &tagged_item, STACK_ITEM_SF);
            
            SFnc *sf = &tagged_item.item.sf;
            SFnc sfd3 = desymmetrize_schofield(*sf);
            free_sfnc(*sf);

            stack_push_with_type(&stack, (void*) &sfd3, STACK_ITEM_SF, pc_loc); 
        
        } else if (strcasecmp(funcname, "WRITE_CF") == 0) {
            Tagged_Stack_Item tagged_item = stack_pop_with_type(&stack, *pc_loc);
            expect_item_on_stack(pc_loc, &tagged_item, STACK_ITEM_CF); 
            
            expect_arg(func);
            const char *filename = processing_params->fs.items[pc].arg;
            CFnc *cf = &tagged_item.item.cf; 

            if (!write_correlation_function(filename, *cf)) {
                PRINT0("ERROR: could not write to file '%s'\n", filename);
                exit(1);
            }

            free_cfnc(*cf);
         
        } else if (strcasecmp(funcname, "WRITE_SF") == 0) {
            Tagged_Stack_Item tagged_item = stack_pop_with_type(&stack, *pc_loc);
            expect_item_on_stack(pc_loc, &tagged_item, STACK_ITEM_SF); 

            expect_arg(func);
            const char *filename = processing_params->fs.items[pc].arg;
            SFnc *sf = &tagged_item.item.sf; 
            
            if (!write_spectral_function(filename, *sf)) {
                PRINT0("ERROR: could not write to file '%s'\n", filename);
                exit(1);
            }

           free_sfnc(*sf);

        } else if (strcasecmp(funcname, "WRITE_SPECTRUM") == 0) { 
            Tagged_Stack_Item tagged_item = stack_pop_with_type(&stack, *pc_loc);
            expect_item_on_stack(pc_loc, &tagged_item, STACK_ITEM_SPECTRUM); 
            
            expect_arg(func);
            const char *filename = processing_params->fs.items[pc].arg;
            Spectrum *sp = &tagged_item.item.sp;

            if (processing_params->spectrum_frequency_max > 0) {
                assert(sp->len > 1);
                double dnu = sp->nu[1] - sp->nu[0];
                size_t n = (size_t) (processing_params->spectrum_frequency_max / dnu);
                sp->len = (sp->len > n) ? n : sp->len;
                
                INFO("Truncating spectrum at max frequency = %.4e cm-1 (resulting npoints: %zu)\n", n*dnu, n);
            }

            if (!write_spectrum(filename, *sp)) {
                PRINT0("ERROR: could not write to file '%s'\n", filename);
                exit(1);
            }

            free_spectrum(*sp); 

        } else {
            PRINT0("\n\n");
            PRINT0("----------------------------------------\n");
            PRINT0("ERROR: running '%s' is not implemented\n", processing_params->fs.items[pc].name);
            PRINT0("----------------------------------------\n");
            assert(false);
        }
        
        _print0_margin = 0;
    }

    if (stack.count > 0) {
        PRINT0("\n\n");
        PRINT0("WARNING: Stack is not empty at the end of processing.\n");
        PRINT0("  Stack trace:\n");
        for (size_t i = 0; i < stack.count; ++i) {
            Loc *loc = &stack.items[i].loc;
            PRINT0("    %zu: %s created at %s:%d:%d\n", 
                   i, STACK_ITEM_TYPES[stack.items[i].typ], loc->input_path, loc->line_number, loc->line_offset);

            switch (stack.items[i].typ) {
                case STACK_ITEM_CF:       free_cfnc(stack.items[i].item.cf); break;
                case STACK_ITEM_SF:       free_sfnc(stack.items[i].item.sf); break;
                case STACK_ITEM_SPECTRUM: free_spectrum(stack.items[i].item.sp); break;
            }
        } 
    } 

    return true;
    
    //if (processing_params->cf_filename == NULL) {
    //    PRINT0("ERROR: required field missing: '%s'\n", KEYWORDS[KEYWORD_CF_FILENAME]);
    //    exit(1);
    //} 
}

void *load_symbol(void *handle, const char *symbol_name, bool allow_undefined) 
{
    void *symbol = dlsym(handle, symbol_name);

    char *err = dlerror();
    if (!allow_undefined && err != NULL) {
        PRINT0("ERROR: could not load '%s': %s\n", symbol_name, err);
        dlclose(handle);
        exit(1);
    }

    return symbol;
}

void setup_dipole(const char *filepath, dipolePtr *dipole_func)
{
    PRINT0("\n\n");
    PRINT0("*****************************************************\n");
    PRINT0("Loading dipole from shared library: %s\n", filepath);

    void *so_handle = dlopen(filepath, RTLD_LAZY);
    if (!so_handle) {
        PRINT0("ERROR: dlopen: %s\n", dlerror());
        exit(1);
    }

    dlerror(); // clear error log
               
    bool allow_undefined = true;
    bool must_be_defined = false;

    void (*dipole_init)(bool);
    dipole_init = load_symbol(so_handle, "dipole_init", allow_undefined);

    if (dipole_init != NULL) {
        if (_wrank == 0) {
            bool log = true;
            printf("INFO: found init function for dipole. Initializing...\n");
            dipole_init(log); 
        } else {
            bool log = false;
            dipole_init(log);
        }
    } else {
        PRINT0("INFO: no 'dipole_init' function found\n");
    }
    
    *dipole_func = (dipolePtr) load_symbol(so_handle, "dipole_lab", must_be_defined);

    PRINT0("Successfully loaded\n");
    PRINT0("*****************************************************\n");
    PRINT0("\n");
}

void setup_pes(InputBlock *input_block) 
{
    PRINT0("\n");
    PRINT0("*****************************************************\n");
    PRINT0("Loading potential from shared library: %s\n", input_block->so_potential);

    void *so_handle = dlopen(input_block->so_potential, RTLD_LAZY);
    if (!so_handle) {
        PRINT0("ERROR: dlopen: %s\n", dlerror());
        exit(1);
    }

    dlerror();
    
    bool allow_undefined = true;
    
    void (*pes_init)(void);
    pes_init = load_symbol(so_handle, "pes_init", allow_undefined);
    
    if (pes_init != NULL) {
        if (_wrank == 0) {
            printf("INFO: found init function for pes. Initializing...\n");
            pes_init(); 
        } else {
            pes_init();
        }
    } else {
        PRINT0("INFO: no 'pes_init' function found\n");
    }
    

    bool must_be_defined = false;
    pes = (pesPtr) load_symbol(so_handle, "pes_lab", must_be_defined);
    dpes = (dpesPtr) load_symbol(so_handle, "dpes_lab", must_be_defined);
    
    PRINT0("Successfully loaded\n");
    PRINT0("*****************************************************\n");
    PRINT0("\n\n");
}

char* shift(int *argc, char ***argv)
{
    assert(*argc > 0);
    char *result = *argv[0];

    *argc -= 1;
    *argv += 1;

    return result; 
}

void usage(const char *program_path)
{
    PRINT0("Usage: %s <configuration-file>\n", program_path);
}

int main(int argc, char* argv[]) 
{
    MPI_Init(&argc, &argv);
    
    INIT_WRANK;
    INIT_WSIZE;

    char *program_path = shift(&argc, &argv);

    if (argc <= 0) {
        usage(program_path);
        PRINT0("ERROR: a configuration file must be provided\n");
        exit(1); 
    }

    const char *filename = shift(&argc, &argv);

    String_Builder sb = {0};
    if (!read_entire_file(filename, &sb)) {
        exit(1);
    }

    PRINT0("*****************************************************\n");
    PRINT0("*********** HAWAII HYBRID v0.1 **********************\n");
    PRINT0("*****************************************************\n");
    PRINT0("Loaded configuration file: %s\n", filename);
    PRINT0("%s\n", sb.items);

    Lexer l = lexer_new(filename, sb.items, sb.items + sb.count);
    //print_lexemes(&l);

    InputBlock input_block = {0}; 
    Monomer monomer1  = {0};
    Monomer monomer2  = {0};
    CalcParams calc_params = {0};
    Processing_Params processing_params = {0};
    parse_params(&l, &calc_params, &input_block, &monomer1, &monomer2, &processing_params);
    
    /*
    print_input_block(&input_block);
    print_params(&params); 
    print_input_block(&input_block);
    print_monomer(&monomer1);
    print_monomer(&monomer2);
    */

    switch (calc_params.calculation_type) {
        case CALCULATION_PR_MU: {
            if (strcmp(input_block.so_dipole_1, input_block.so_dipole_2) != 0) {
                PRINT0("ERROR: for PR_MU calculation we expect to have only one dipole\n");
                exit(1);
            }
            
            assert(input_block.so_dipole_1 != NULL);
            setup_dipole(input_block.so_dipole_1, &dipole_1);
            setup_pes(&input_block);

            MoleculeSystem ms = init_ms_from_monomers(input_block.reduced_mass, &monomer1, &monomer2, 0);
            SFnc sf = calculate_spectral_function_using_prmu_representation_and_save(&ms, &calc_params, input_block.Temperature);

            UNUSED(sf);
            break; 
        }
        case CALCULATION_CORRELATION_SINGLE: {
            setup_dipole(input_block.so_dipole_1, &dipole_1);
            if (input_block.so_dipole_2 != NULL) {
                setup_dipole(input_block.so_dipole_2, &dipole_2);
            } else {
                dipole_2 = dipole_1;
            }

            setup_pes(&input_block);

            MoleculeSystem ms = init_ms_from_monomers(input_block.reduced_mass, &monomer1, &monomer2, 0);
            calculate_correlation_and_save(&ms, &calc_params, input_block.Temperature);
           
            if (_wrank == 0) { 
                if (!run_processing(&processing_params)) {
                    PRINT0("ERROR: could not run commands in PROCESSING block\n");
                    exit(1); 
                }
            }

            break;
        } 
        case CALCULATION_CORRELATION_ARRAY: {
            setup_dipole(input_block.so_dipole_1, &dipole_1);
            
            if (input_block.so_dipole_2 != NULL) {
                setup_dipole(input_block.so_dipole_2, &dipole_2);
            } else {
                dipole_2 = dipole_1;
            }

            setup_pes(&input_block);

            MoleculeSystem ms = init_ms_from_monomers(input_block.reduced_mass, &monomer1, &monomer2, 0);

            double base_temperature = input_block.Temperature;
            CFncArray ca = calculate_correlation_array_and_save(&ms, &calc_params, base_temperature);

            UNUSED(ca);
            break;
        }
        case CALCULATION_PROCESSING: {
            if ((_wrank == 0) && !run_processing(&processing_params)) {
                PRINT0("ERROR: could not run commands in PROCESSING block\n");
                exit(1); 
            }

            break; 
        }
        case CALCULATION_PHASE_SPACE_M2: {
            setup_dipole(input_block.so_dipole_1, &dipole_1);
            if (input_block.so_dipole_2 != NULL) {
                setup_dipole(input_block.so_dipole_2, &dipole_2);
            } else {
                dipole_2 = dipole_1;
            }
            
            setup_pes(&input_block);
            
            MoleculeSystem ms = init_ms_from_monomers(input_block.reduced_mass, &monomer1, &monomer2, 0);
            
            size_t hep_m2_niterations = 12;
            if (calc_params.hep_m2_niterations > 0) hep_m2_niterations = calc_params.hep_m2_niterations;
            
            size_t hep_m2_npoints = 1000000;
            if (calc_params.hep_m2_npoints > 0) hep_m2_npoints = calc_params.hep_m2_npoints;
            
            if (calc_params.sampler_Rmin <= 0.0) {
                PRINT0("ERROR: minimum value for R throught integration (sampler_Rmin) must be set\n");
                exit(1);
            }
            
            if (calc_params.sampler_Rmax <= 0.0) {
                PRINT0("ERROR: maximum value for R throught integration (sampler_Rmin) must be set\n");
                exit(1);
            }

            if (input_block.Temperature > 0) {
                double T = input_block.Temperature; 
                double pf_analytic = analytic_full_partition_function_by_V(&ms, T);
                INFO("Analytic partition function divided by V: %.5e\n", pf_analytic);

                double hep_M2, hep_M2_err; 
                c_mpi_perform_integration(&ms, INTEGRAND_M2, &calc_params, T, hep_m2_niterations, hep_m2_npoints, &hep_M2, &hep_M2_err);

                hep_M2     *= SecondCoeff / pf_analytic;
                hep_M2_err *= SecondCoeff / pf_analytic;
                INFO("T = %.2e => M2: %.5e\n", T, hep_M2);
            } else {
                PRINT0("ERROR: no temperature is provided to run CALCULATION_PHASE_SPACE_M2\n");
                exit(1); 
            }

            break;
        }
        case CALCULATION_PHASE_SPACE_M0: {
            setup_dipole(input_block.so_dipole_1, &dipole_1);
            if (input_block.so_dipole_2 != NULL) {
                setup_dipole(input_block.so_dipole_2, &dipole_2);
            } else {
                dipole_2 = dipole_1;
            }

            setup_pes(&input_block);

            MoleculeSystem ms = init_ms_from_monomers(input_block.reduced_mass, &monomer1, &monomer2, 0);
        
            size_t hep_m0_niterations = 12;
            if (calc_params.hep_m0_niterations > 0) hep_m0_niterations = calc_params.hep_m0_niterations;

            size_t hep_m0_npoints = 1000000;
            if (calc_params.hep_m0_npoints > 0) hep_m0_npoints = calc_params.hep_m0_npoints;

            if (calc_params.sampler_Rmin <= 0.0) {
                PRINT0("ERROR: minimum value for R throught integration (sampler_Rmin) must be set\n");
                exit(1);
            }
            
            if (calc_params.sampler_Rmax <= 0.0) {
                PRINT0("ERROR: maximum value for R throught integration (sampler_Rmin) must be set\n");
                exit(1);
            }

            if (input_block.num_temperatures > 0) {
                for (size_t i = 0; i < input_block.num_temperatures; ++i) {
                    double T = input_block.temperatures[i]; 
                    double pf_analytic = analytic_full_partition_function_by_V(&ms, T);
                    INFO("Analytic partition function divided by V: %.5e\n", pf_analytic);
                    
                    double hep_M0, hep_M0_err; 
                    c_mpi_perform_integration(&ms, INTEGRAND_M0, &calc_params, T, hep_m0_niterations, hep_m0_npoints, &hep_M0, &hep_M0_err);
                    
                    hep_M0     *= ZeroCoeff / pf_analytic;
                    hep_M0_err *= ZeroCoeff / pf_analytic;
                    INFO("T = %.2e => M0: %.5e\n", T, hep_M0);
                }
            } else if (input_block.Temperature > 0) {
                double T = input_block.Temperature; 
                double pf_analytic = analytic_full_partition_function_by_V(&ms, T);
                INFO("Analytic partition function divided by V: %.5e\n", pf_analytic);

                double hep_M0, hep_M0_err; 
                c_mpi_perform_integration(&ms, INTEGRAND_M0, &calc_params, T, hep_m0_niterations, hep_m0_npoints, &hep_M0, &hep_M0_err);

                hep_M0     *= ZeroCoeff / pf_analytic;
                hep_M0_err *= ZeroCoeff / pf_analytic;
                INFO("T = %.2e => M0: %.5e\n", T, hep_M0);
            } else {
                PRINT0("ERROR: no temperature is provided to run CALCULATION_PHASE_SPACE_M0\n");
                exit(1); 
            }

            break;
        }
        case CALCULATION_NONE: UNREACHABLE(""); 
        case CALCULATION_TYPES_COUNT: UNREACHABLE(""); 
    } 

    /*
    double *q = malloc(ms.Q_SIZE*sizeof(double));
    q[2] = 6.0;
    q[0] = 0.1;
    q[1] = 0.2;
    q[3] = 0.3;
    q[4] = 0.4;
    double dip[3];
    dipole(q, dip);  
    printf("DIPOLE VALUE: %.5e %.5e %.5e\n", dip[0], dip[1], dip[2]);
    */
    
    MPI_Finalize();

    return 0; 
}
