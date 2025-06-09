#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <inttypes.h>

#include "hawaii.h"

#define BUFF_SIZE 256

// TODO: define those as function pointers in 'hawaii.h'
double pes(double *q) {
    UNUSED(q);
    assert(false);
}

void dpes(double *q, double *dVdq) {
    UNUSED(q);
    UNUSED(dVdq);
    assert(false);
}

const char *PAIR_STATES[] = {
    "FREE_AND_METASTABLE",
    "BOUND",
};


const char* CALCULATION_TYPES[] = {
    "PR_MU",
    "CORRELTION_SINGLE",
    "CORRELATION_ARRAY",
};

#define BLOCKS_MAX 16
typedef enum {
    BLOCK_NONE,
    BLOCK_INPUT,
    BLOCK_MONOMER,
} Block_Name;

typedef struct {
    Block_Name items[BLOCKS_MAX];
    size_t count;
} Block_Names; 

static_assert(MONOMER_COUNT == 6);
MonomerType MONOMER_TYPES[MONOMER_COUNT] = {
    ATOM, LINEAR_MOLECULE, LINEAR_MOLECULE_REQ_INTEGER, LINEAR_MOLECULE_REQ_HALFINTEGER, ROTOR, ROTOR_REQUANTIZED_ROTATION,
};

const char* display_monomer_type(MonomerType t) {
    switch (t) {
        case ATOM:                            return "ATOM";
        case LINEAR_MOLECULE:                 return "LINEAR_MOLECULE";
        case LINEAR_MOLECULE_REQ_INTEGER:     return "LINEAR_MOLECULE_REQ_INTEGER";
        case LINEAR_MOLECULE_REQ_HALFINTEGER: return "LINEAR_MOLECULE_REQ_HALFINTEGER";
        case ROTOR:                           return "ROTOR";
        case ROTOR_REQUANTIZED_ROTATION:      return "ROTOR_REQUANTIZED_ROTATION";
    }

    return NULL;
}

typedef struct {
    MonomerType t;
    double II[3];
    double DJ;
} MonomerBlock;

typedef struct {
    double reduced_mass;
    const char* so_potential;
    const char* so_dipole;
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
    TOKEN_OCURLY,
    TOKEN_CCURLY,
    TOKEN_COMMA,
    TOKEN_EQ,
    TOKEN_COUNT, // sentinel value for array size
} Token_Type;

const char *TOKEN_TYPES[TOKEN_COUNT] = {
      [TOKEN_EOF]      = "EOF",
      [TOKEN_BLOCK]    = "BLOCK",
      [TOKEN_KEYWORD]  = "KEYWORD",
      [TOKEN_STRING]   = "STRING",
      [TOKEN_DQSTRING] = "DOUBLE-QUOTED STRING",
      [TOKEN_INTEGER]  = "INTEGER",
      [TOKEN_FLOAT]    = "FLOAT",
      [TOKEN_OCURLY]   = "{",
      [TOKEN_CCURLY]   = "}",
      [TOKEN_COMMA]    = ",",
      [TOKEN_EQ]       = "=",
}; 

static_assert(TOKEN_COUNT == 11);

const char *PUNCTS[TOKEN_COUNT] = {
    [TOKEN_COMMA]  = ",",
    [TOKEN_OCURLY] = "{",
    [TOKEN_CCURLY] = "}",
    [TOKEN_EQ]     = "=",
};

typedef struct {
    Token_Type t;
    char *s;
    Loc loc;
} Token; 

typedef enum {
    /* INPUT BLOCK */
    KEYWORD_CALCULATION_TYPE,
    KEYWORD_PAIR_STATE,
    KEYWORD_REDUCED_MASS,
    KEYWORD_SO_POTENTIAL,
    KEYWORD_SO_DIPOLE,
    KEYWORD_NITERATIONS,
    KEYWORD_TOTAL_TRAJECTORIES,
    KEYWORD_CVODE_TOLERANCE,
    KEYWORD_SAMPLING_TIME,
    KEYWORD_MAXTRAJECTORYLENGTH,
    KEYWORD_INITIALM0_NPOINTS,
    KEYWORD_INITIALM2_NPOINTS,
    KEYWORD_R0,
    KEYWORD_APPROXIMATEFREQUENCYMAX,
    KEYWORD_TORQUE_CACHE_LEN,
    KEYWORD_TORQUE_LIMIT,
    KEYWORD_JINI_HISTOGRAM_BINS,
    KEYWORD_JINI_HISTOGRAM_MAX,
    KEYWORD_JFIN_HISTOGRAM_BINS,
    KEYWORD_JFIN_HISTOGRAM_MAX,
    KEYWORD_ORTHO_STATE_WEIGHT,
    KEYWORD_PARA_STATE_WEIGHT,
    /* MONOMER BLOCK */
    KEYWORD_MONOMER_TYPE,
    KEYWORD_DJ,
    KEYWORD_II,
    KEYWORD_COUNT,
} Keyword;

const char* KEYWORDS[KEYWORD_COUNT] = {
    [KEYWORD_CALCULATION_TYPE]        = "CALCULATION_TYPE",
    [KEYWORD_PAIR_STATE]              = "PAIR_STATE",
    [KEYWORD_REDUCED_MASS]            = "REDUCED_MASS",
    [KEYWORD_SO_POTENTIAL]            = "SO_POTENTIAL",
    [KEYWORD_SO_DIPOLE]               = "SO_DIPOLE",
    [KEYWORD_NITERATIONS]             = "NITERATIONS",
    [KEYWORD_TOTAL_TRAJECTORIES]      = "TOTAL_TRAJECTORIES",
    [KEYWORD_CVODE_TOLERANCE]         = "CVODE_TOLERANCE",
    [KEYWORD_SAMPLING_TIME]           = "SAMPLING_TIME",
    [KEYWORD_MAXTRAJECTORYLENGTH]     = "MAXTRAJECTORYLENGTH",
    [KEYWORD_INITIALM0_NPOINTS]       = "INITIALM0_NPOINTS",
    [KEYWORD_INITIALM2_NPOINTS]       = "INITIALM2_NPOINTS",
    [KEYWORD_R0]                      = "R0",
    [KEYWORD_APPROXIMATEFREQUENCYMAX] = "APPROXIMATEFREQUENCYMAX",
    [KEYWORD_TORQUE_CACHE_LEN]        = "TORQUE_CACHE_LEN",
    [KEYWORD_TORQUE_LIMIT]            = "TORQUE_LIMIT",
    [KEYWORD_JINI_HISTOGRAM_BINS]     = "JINI_HISTOGRAM_BINS",
    [KEYWORD_JINI_HISTOGRAM_MAX]      = "JINI_HISTOGRAM_MAX",
    [KEYWORD_JFIN_HISTOGRAM_BINS]     = "JFIN_HISTOGRAM_BINS",
    [KEYWORD_JFIN_HISTOGRAM_MAX]      = "JFIN_HISTOGRAM_MAX",
    [KEYWORD_ORTHO_STATE_WEIGHT]      = "ORTHO_STATE_WEIGHT",
    [KEYWORD_PARA_STATE_WEIGHT]       = "PARA_STATE_WEIGHT",
    [KEYWORD_MONOMER_TYPE]            = "MONOMER_TYPE",
    [KEYWORD_DJ]                      = "DJ",
    [KEYWORD_II]                      = "II",
}; 
static_assert(KEYWORD_COUNT == 25);

Token_Type EXPECT_TOKEN[KEYWORD_COUNT] = {
    [KEYWORD_CALCULATION_TYPE]        = TOKEN_STRING,
    [KEYWORD_PAIR_STATE]              = TOKEN_STRING,
    [KEYWORD_REDUCED_MASS]            = TOKEN_FLOAT,
    [KEYWORD_SO_POTENTIAL]            = TOKEN_DQSTRING,
    [KEYWORD_SO_DIPOLE]               = TOKEN_DQSTRING,
    [KEYWORD_NITERATIONS]             = TOKEN_INTEGER,
    [KEYWORD_TOTAL_TRAJECTORIES]      = TOKEN_INTEGER,
    [KEYWORD_CVODE_TOLERANCE]         = TOKEN_FLOAT,
    [KEYWORD_SAMPLING_TIME]           = TOKEN_FLOAT,
    [KEYWORD_MAXTRAJECTORYLENGTH]     = TOKEN_INTEGER,
    [KEYWORD_INITIALM0_NPOINTS]       = TOKEN_INTEGER,
    [KEYWORD_INITIALM2_NPOINTS]       = TOKEN_INTEGER,
    [KEYWORD_R0]                      = TOKEN_FLOAT,
    [KEYWORD_APPROXIMATEFREQUENCYMAX] = TOKEN_FLOAT,
    [KEYWORD_TORQUE_CACHE_LEN]        = TOKEN_INTEGER,
    [KEYWORD_TORQUE_LIMIT]            = TOKEN_FLOAT,
    [KEYWORD_JINI_HISTOGRAM_BINS]     = TOKEN_INTEGER,
    [KEYWORD_JINI_HISTOGRAM_MAX]      = TOKEN_FLOAT,
    [KEYWORD_JFIN_HISTOGRAM_BINS]     = TOKEN_INTEGER,
    [KEYWORD_JFIN_HISTOGRAM_MAX]      = TOKEN_FLOAT,
    [KEYWORD_ORTHO_STATE_WEIGHT]      = TOKEN_FLOAT,
    [KEYWORD_PARA_STATE_WEIGHT]       = TOKEN_FLOAT,
    [KEYWORD_MONOMER_TYPE]            = TOKEN_STRING,
    [KEYWORD_DJ]                      = TOKEN_FLOAT,
    [KEYWORD_II]                      = TOKEN_OCURLY,
};

/*
int trim_whitespace(char *s)
{
    if (!s) return 0;
   
    int trimmed_from_beginning = 0; 
    char *start = s;
    while (isspace(*start)) { 
        start++;
        trimmed_from_beginning++;
    }

    char *end = s + strlen(s) - 1;
    while (end > start && isspace(*end)) end--;

    *(end + 1) = '\0';
    if (start != s) memmove(s, start, end - start + 2); 

    return trimmed_from_beginning; 
}

bool parse_key_value(char *line, Token *key, Token *value)
{
    char *equals = strchr(line, '=');
    if (!equals) {
        printf("ERROR: %s:%d: expected '=' but not found\n", key->loc.input_path, key->loc.line_number);
        return false;
    }
    
    (*value).s = equals + 1;
    (*key).s = line;
    *equals = '\0';
  
    trim_whitespace((*key).s);
    value->loc.line_offset += (int) ((*value).s - line);
    value->loc.line_offset += trim_whitespace((*value).s);

    return true;
}

int expect_int(char *value, Loc loc)
{
    char *endptr;
    long p = strtol(value, &endptr, 10);
    if (endptr != value + strlen(value)) {
        printf("ERROR: %s:%d:%d: could not parse string '%s' as integer\n", loc.input_path, loc.line_number, loc.line_offset, value);
        exit(1); 
    }

    return (int) p;
}

double expect_double(char* value, Loc *loc)
{
    char *endptr;
    double p = strtof(value, &endptr);
    if (endptr != value + strlen(value)) {
        printf("ERROR: %s:%d:%d: could not parse string '%s' as float\n", loc->input_path, loc->line_number, loc->line_offset, value);
        exit(1);
    }
   
    loc->line_offset += strlen(value);

    return p;    
}

double *expect_double_array(char *value, Loc loc)
{
    if (value[0] != '{') { 
        printf("ERROR: %s:%d:%d: expected opening curly-brace to start the array\n", loc.input_path, loc.line_number, loc.line_offset);
        exit(1); 
    }

    // advance over '{' 
    value = value + 1; 
    loc.line_offset++; 

    double *a = (double*) malloc(3*sizeof(double));
    memset(a, 0, 3*sizeof(double));

    size_t c = 0;

    do {
        int n = trim_whitespace(value);
        loc.line_offset += n;

        a[c++] = expect_double(value, &loc);
        loc.line_offset += n;

        assert(false);
    } while(true);

    if (value[strlen(value) - 1] != '}') {
        printf("ERROR: %s:%d:%d: expected closing curly-brace at the end of the array\n", loc.input_path, loc.line_number, loc.line_offset);
    }

    return a;
}

char *expect_dqstring(char *value, Loc loc) 
{
    if (value[0] != '\"' || value[strlen(value) - 1] != '\"') {
        printf("ERROR: %s:%d:%d: expected a double-quoted string\n", loc.input_path, loc.line_number, loc.line_offset);
        exit(1);
    }

    return strdup(value);
}

char *expect_string(char *value, Loc loc)
{
    UNUSED(loc);
    return strdup(value);
}
                
void parse_monomer_block_line(char *line, MonomerBlock *monomer, Loc loc)
{
    Token key = {
        .s = NULL,
        .loc = loc,
    };

    Token value = {
        .s = NULL,
        .loc = loc,
    };
    
    if (!parse_key_value(line, &key, &value)) return;
    
    printf("key = '%s', value = '%s'\n", key.s, value.s);

    if (strcasecmp(key.s, "MONOMER_TYPE") == 0) {
        char *monomer_type_str = expect_string(value.s, value.loc);

        if (strcasecmp(monomer_type_str, "ATOM") == 0) {
            monomer->t = ATOM;
        } else if (strcasecmp(monomer_type_str, "LINEAR_MOLECULE") == 0) {
            monomer->t = LINEAR_MOLECULE;
        } else if (strcasecmp(monomer_type_str, "LINEAR_MOLECULE_REQ_INTEGER") == 0) {
            monomer->t = LINEAR_MOLECULE_REQ_INTEGER;
        } else if (strcasecmp(monomer_type_str, "LINEAR_MOLECULE_REQ_HALFINTEGER") == 0) {
            monomer->t = LINEAR_MOLECULE_REQ_HALFINTEGER;
        } else if (strcasecmp(monomer_type_str, "ROTOR") == 0) {
            monomer->t = ROTOR;
        } else {
            printf("ERROR: %s:%d:%d: unknown monomer type '%s'\n", value.loc.input_path, value.loc.line_number, value.loc.line_offset, monomer_type_str);
            printf("Available monomer types:\n");
            for (size_t i = 0; i < sizeof(MONOMER_TYPES)/sizeof(MONOMER_TYPES[0]); ++i) {
                printf(" %s\n", MONOMER_TYPES[i]); 
            }
            exit(1);
        }
    } else if (strcasecmp(key.s, "DJ") == 0) {
        monomer->DJ = expect_double(value.s, &value.loc);
    } else if (strcasecmp(key.s, "II") == 0) {
        double *a = expect_double_array(value.s, value.loc);
        memcpy(monomer->II, a, 3*sizeof(double));
        free(a); 
    } else {
        printf("ERROR: %s:%d:%d: unknown key '%s'\n", key.loc.input_path, key.loc.line_number, key.loc.line_offset, key.s);
        exit(1);
    }
}

void parse_input_block_line(char *line, InputBlock *input_block, CalcParams *params, Loc loc)
{
    Token key = {
        .s = NULL,
        .loc = loc,
    };

    Token value = {
        .s = NULL,
        .loc = loc,
    };

    // printf("offset = %d\n", key.loc.line_offset);
    if (!parse_key_value(line, &key, &value)) return;

    //printf("key = '%s', value = '%s'\n", key.s, value.s);

    if (strcasecmp(key.s, "CALCULATION_TYPE") == 0) {
        char *calculation_type_str = expect_string(value.s, value.loc);

        if (strcasecmp(calculation_type_str, "PR_MU") == 0) {
            input_block->calculation_type = PR_MU;
        } else if (strcasecmp(calculation_type_str, "CORRELATION_SINGLE") == 0) {
            input_block->calculation_type = CORRELATION_SINGLE;
        } else if (strcasecmp(calculation_type_str, "CORRELATION_ARRAY") == 0) {
            input_block->calculation_type = CORRELATION_ARRAY;
        } else {
            printf("ERROR: %s:%d:%d: unknown calculation type '%s'\n", value.loc.input_path, value.loc.line_number, value.loc.line_offset, calculation_type_str);
            printf("Available calculation types:\n");
            for (size_t i = 0; i < sizeof(CALCULATION_TYPES)/sizeof(CALCULATION_TYPES[0]); ++i) {
                printf("  %s\n", CALCULATION_TYPES[i]);
            }
            exit(1);
        }
        
        free(calculation_type_str);

    } else if (strcasecmp(key.s, "REDUCED_MASS") == 0) {
        input_block->reduced_mass = expect_double(value.s, &value.loc);
    } else if (strcasecmp(key.s, "SO_POTENTIAL") == 0) {
        input_block->so_potential = expect_dqstring(value.s, value.loc); 
    } else if (strcasecmp(key.s, "SO_DIPOLE") == 0) {
        input_block->so_dipole = expect_dqstring(value.s, value.loc);
    } else if (strcasecmp(key.s, "II1") == 0) {
    } else if (strcasecmp(key.s, "PAIR_STATE") == 0) {
        char *pair_state_str = expect_string(value.s, value.loc);
        
        if (strcasecmp(pair_state_str, "FREE_AND_METASTABLE") == 0) {
            params->ps = FREE_AND_METASTABLE;
        } else if (strcasecmp(pair_state_str, "BOUND") == 0) {
            params->ps = BOUND;
        }

        free(pair_state_str);

    } else if (strcasecmp(key.s, "NITERATIONS") == 0) {
        params->niterations = expect_int(value.s, value.loc); 
    } else if (strcasecmp(key.s, "TOTAL_TRAJECTORIES") == 0) {
        params->total_trajectories = expect_int(value.s, value.loc); 
    } else if (strcasecmp(key.s, "CVODE_TOLERANCE") == 0) {
        params->cvode_tolerance = expect_double(value.s, &value.loc); 
    } else if (strcasecmp(key.s, "SAMPLING_TIME") == 0) {
        params->sampling_time = expect_double(value.s, &value.loc); 
    } else if (strcasecmp(key.s, "MAXTRAJECTORYLENGTH") == 0) {
        params->MaxTrajectoryLength = expect_int(value.s, value.loc); 
    } else if (strcasecmp(key.s, "INITIALM0_NPOINTS") == 0) {
        params->initialM0_npoints = expect_int(value.s, value.loc); 
    } else if (strcasecmp(key.s, "INITIALM2_NPOINTS") == 0) {
        params->initialM2_npoints = expect_int(value.s, value.loc); 
    } else if (strcasecmp(key.s, "R0") == 0) {
        params->R0 = expect_double(value.s, &value.loc); 
    } else if (strcasecmp(key.s, "APPROXIMATEFREQUENCYMAX") == 0) {
        params->ApproximateFrequencyMax = expect_double(value.s, &value.loc); 
    } else if (strcasecmp(key.s, "TORQUE_CACHE_LEN") == 0) {
        params->torque_cache_len = expect_int(value.s, value.loc); 
    } else if (strcasecmp(key.s, "TORQUE_LIMIT") == 0) {
        params->torque_limit = expect_double(value.s, &value.loc); 
    } else if (strcasecmp(key.s, "JINI_HISTOGRAM_BINS") == 0) {
        params->jini_histogram_bins = expect_int(value.s, value.loc); 
    } else if (strcasecmp(key.s, "JINI_HISTOGRAM_MAX") == 0) {
        params->jini_histogram_max = expect_double(value.s, &value.loc); 
    } else if (strcasecmp(key.s, "JFIN_HISTOGRAM_BINS") == 0) {
        params->jfin_histogram_bins = expect_int(value.s, value.loc); 
    } else if (strcasecmp(key.s, "JFIN_HISTOGRAM_MAX") == 0) {
        params->jfin_histogram_max = expect_double(value.s, &value.loc); 
    } else if (strcasecmp(key.s, "ORTHO_STATE_WEIGHT") == 0) {
        params->ortho_state_weight = expect_double(value.s, &value.loc); 
    } else if (strcasecmp(key.s, "PARA_STATE_WEIGHT") == 0) {
        params->para_state_weight = expect_double(value.s, &value.loc); 
    } else {
        printf("ERROR: %s:%d:%d: unknown key '%s'\n", key.loc.input_path, key.loc.line_number, key.loc.line_offset, key.s);
        exit(1);
    } 
}
*/

typedef struct {
    char *input_stream;
    char *parse_point;
    char *eof;
    Token_Type token_type;
    String_Builder string_storage;
    int64_t int_number;
    double double_number;
    Keyword keyword_type;
    Loc loc;
} Lexer;

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
    return isalnum(c) || c == '_';
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

    char c = peek_char(l);;
    if (c == '\0') {
        l->token_type = TOKEN_EOF;
        return false; 
    }

    for (size_t i = 0; i < sizeof(PUNCTS)/sizeof(PUNCTS[0]); ++i) {
        const char *prefix = PUNCTS[i];
        if (skip_prefix(l, prefix)) {
            l->token_type = (Token_Type) i; 
        } 
    } 

    if (c == '&') { // block begin
        (*l).token_type = TOKEN_BLOCK;
        (*l).string_storage.count = 0;
        skip_char(l); // consume '&'

        for ( c = peek_char(l); c != '\0'; c = peek_char(l)) {
            if (is_identifier(c)) {
                da_append(&(*l).string_storage, c);
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
        
        for (c = peek_char(l); is_identifier(c); c = peek_char(l)) {
            da_append(&(*l).string_storage, c);
            skip_char(l);
        }

        da_append(&(*l).string_storage, 0);

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
            
                l->double_number += fraction/divisor;
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
            printf("ERROR: %s:%d:%d: unfinished string literal\n", l->loc.input_path, l->loc.line_number, l->loc.line_offset);
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

/*
void parse_params_from_file(const char *filename, InputBlock *input_block, MonomerBlock *monomer1_block, MonomerBlock *monomer2_block, CalcParams *params)
{
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("ERROR: could not open the file '%s'\n", filename);
        exit(1);
    }

    char line[BUFF_SIZE];

    Loc loc = {
        .input_path = filename,
        .line_number = 0,
    };

    Block_Name current_block = BLOCK_NONE;
    Block_Names blocks_found = {0};
        
    while (fgets(line, sizeof(line), fp)) 
    {
        char *ptr = line;
        loc.line_offset = 1;
        loc.line_number++;

        while (isspace((unsigned char) *ptr)) {
            loc.line_offset++; 
            ptr++;
        }
        
        if ((*ptr == '!') || (*ptr == '\0')) continue; 
      

        if (strncmp(ptr, "&INPUT", strlen("&INPUT")) == 0) {
            if (current_block != BLOCK_NONE) {
                printf("ERROR: %s:%d:%d: blocks cannot be nested\n", loc.input_path, loc.line_number, loc.line_offset);
                exit(1);
            }

            for (size_t i = 0; i < blocks_found.count; ++i) {
                if (blocks_found.items[i] == BLOCK_INPUT) {
                    printf("ERROR: %s:%d:%d: $INPUT block cannot be defined multiple times\n", loc.input_path, loc.line_number, loc.line_offset);
                    exit(1);
                } 
            }

            current_block = BLOCK_INPUT;
            blocks_found.items[blocks_found.count++] = BLOCK_INPUT;
            continue; 
        }

        if (strncmp(ptr, "&MONOMER", strlen("&MONOMER")) == 0) {
            if (current_block != BLOCK_NONE) {
                printf("ERROR: %s:%d:%d: blocks cannot be nested\n", loc.input_path, loc.line_number, loc.line_offset);
                exit(1);
            }

            size_t monomer_blocks_count = 0;
            for (size_t i = 0; i < blocks_found.count; ++i) {
                if (blocks_found.items[i] == BLOCK_MONOMER) monomer_blocks_count++; 
            }

            if (monomer_blocks_count >= 2) {
                printf("ERROR: %s:%d:%d: only two &MONOMER blocks could be defined\n", loc.input_path, loc.line_number, loc.line_offset);
                exit(1);
            }

            current_block = BLOCK_MONOMER;
            blocks_found.items[blocks_found.count++] = BLOCK_MONOMER;
            continue;
        }

        if (strncmp(ptr, "&END", strlen("&END")) == 0) {
            if (current_block != BLOCK_NONE) {
                current_block = BLOCK_NONE;
                continue;
            } else {
                printf("ERROR: %s:%d:%d: block end marker (&END) encountered without a corresponding block start\n", loc.input_path, loc.line_number, loc.line_offset);
                exit(1); 
            }
        } 

        switch (current_block) { 
            case BLOCK_INPUT: {
                parse_input_block_line(ptr, input_block, params, loc); 
                break;
            }
            case BLOCK_MONOMER: { 
                parse_monomer_block_line(ptr, monomer1_block, loc);
                break;
            } 
        }
    }


    fclose(fp);
}
*/

void print_input_block(InputBlock *input_block) {
    printf("Input Block:\n");
    printf("  reduced_mass = %.5e\n", input_block->reduced_mass);
    printf("  so_potential = %s\n", input_block->so_potential);
    printf("  so_dipole    = %s\n", input_block->so_dipole);
}

void print_monomer(MonomerBlock *monomer_block) {
    printf("Monomer Block:\n");
    printf("  t = %s\n", display_monomer_type(monomer_block->t));
    printf("  I = {%.5e, %.5e, %.5e}\n", monomer_block->II[0], monomer_block->II[1], monomer_block->II[2]);
    printf("  DJ = %.5e\n", monomer_block->DJ);  
}

void print_params(CalcParams *params) {
    printf("Params:\n");
    printf("  pair state              = %s\n", PAIR_STATES[params->ps]);
    printf("  niterations             = %zu\n", params->niterations);
    printf("  total_trajectories      = %zu\n", params->total_trajectories); 
    printf("  cvode_tolerance         = %.5e\n", params->cvode_tolerance);         
    printf("  sampling_time           = %.5e\n", params->sampling_time);
    printf("  MaxTrajectoryLength     = %zu\n", params->MaxTrajectoryLength);
    printf("  initialM0_npoints       = %zu\n", params->initialM0_npoints);
    printf("  initialM2_npoints       = %zu\n", params->initialM2_npoints);
    printf("  R0                      = %.5e\n", params->R0);
    printf("  ApproximateFrequencyMax = %.5e\n", params->ApproximateFrequencyMax);
    printf("  torque_cache_len        = %zu\n", params->torque_cache_len);
    printf("  torque_limit            = %.5e\n", params->torque_limit);
    printf("  jini_histogram_bins     = %zu\n", params->jini_histogram_bins);
    printf("  jini_histogram_max      = %.5e\n", params->jini_histogram_max);
    printf("  jfin_histogram_bins     = %zu\n", params->jfin_histogram_bins);
    printf("  jfin_histogram_max      = %.5e\n", params->jfin_histogram_max);
    printf("  ortho_state_weight      = %.5e\n", params->ortho_state_weight);
    printf("  para_state_weight       = %.5e\n", params->para_state_weight);
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

    sb_reserve(sb, m);
    fread(sb->items + sb->count, m, 1, fp);
    sb->count += m;

    if (ferror(fp)) {
       return_defer(false);
    }

defer:
    if (!result) printf("ERROR: Could not read file %s: %s", path, strerror(errno));
    if (fp) fclose(fp);
    return result;
}

void print_lexeme(Lexer *l) {
    switch (l->token_type) {
      case TOKEN_KEYWORD: {
        printf("KEYWORD: %s\n", l->string_storage.items);
        break;
      }
      case TOKEN_STRING: { 
        printf("%s: '%s'\n", TOKEN_TYPES[l->token_type], l->string_storage.items);
        break;
      }
      case TOKEN_BLOCK: {
        printf("%s: %s\n", TOKEN_TYPES[l->token_type], l->string_storage.items);
        break;
      }
      case TOKEN_INTEGER: {
        printf("%s: %"PRId64"\n", TOKEN_TYPES[l->token_type], l->int_number);
        break;
      }
      case TOKEN_FLOAT: {
        printf("%s: %.16e\n", TOKEN_TYPES[l->token_type], l->double_number);
        break;
      }
      case TOKEN_DQSTRING: {
        printf("%s: \"%s\"\n", TOKEN_TYPES[l->token_type], l->string_storage.items); 
        break;
      }
      case TOKEN_EOF: {
        printf("EOF\n");
        break;
      }
      default: {
        printf("TOKEN: %s\n", TOKEN_TYPES[l->token_type]);
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

void expect_token(Lexer *l, Token_Type expected) {
    Token_Type t = l->token_type;
    print_lexeme(l);

    // TODO: use display_token
    if (t != expected) {
        printf("ERROR: %s:%d:%d: expected '%s' but got type '%s'\n", l->loc.input_path, l->loc.line_number, l->loc.line_offset, 
               TOKEN_TYPES[expected], TOKEN_TYPES[t]);
        exit(1);
    }
}

void get_and_expect_token(Lexer *l, Token_Type token) {
    get_token(l);
    expect_token(l, token);
}

void parse_input_block(Lexer *l, InputBlock *input_block, CalcParams *params) 
{
    while (true) {
        get_token(l);
        if ((l->token_type == TOKEN_BLOCK) && (strcasecmp(l->string_storage.items, "END") == 0)) {
            print_lexeme(l);
            return;
        }

        expect_token(l, TOKEN_KEYWORD);
        Keyword keyword_type = l->keyword_type;
        Token_Type expect_token = EXPECT_TOKEN[keyword_type]; 

        get_and_expect_token(l, TOKEN_EQ);
        get_and_expect_token(l, expect_token);
        
        switch (keyword_type) {
            case KEYWORD_CALCULATION_TYPE: {
                if (strcasecmp(l->string_storage.items, "PR_MU") == 0) {
                    params->calculation_type = PR_MU;
                } else if (strcasecmp(l->string_storage.items, "CORRELATION_SINGLE") == 0) {
                    params->calculation_type = CORRELATION_SINGLE;
                } else if (strcasecmp(l->string_storage.items, "CORRELATION_ARRAY") == 0) {
                    params->calculation_type = CORRELATION_ARRAY;
                } else {
                    printf("ERROR: %s:%d:%d: unknown calculation type '%s'\n", l->loc.input_path, l->loc.line_number, l->loc.line_offset, l->string_storage.items);
                    printf("Available calculation types:\n");
                    
                    for (size_t i = 0; i < sizeof(CALCULATION_TYPES)/sizeof(CALCULATION_TYPES[0]); ++i) {
                        printf("  %s\n", CALCULATION_TYPES[i]);
                    }

                    exit(1);
                }

                break;
            }
            case KEYWORD_PAIR_STATE: {
                if (strcasecmp(l->string_storage.items, "FREE_AND_METASTABLE") == 0) {
                    params->ps = FREE_AND_METASTABLE;
                } else if (strcasecmp(l->string_storage.items, "BOUND") == 0) {
                    params->ps = BOUND;
                } else {
                    printf("ERROR: %s:%d:%d: unknown pair state '%s'\n", l->loc.input_path, l->loc.line_number, l->loc.line_offset, l->string_storage.items);
                    printf("Available pair states:\n");

                    for (size_t i = 0; i < sizeof(PAIR_STATES)/sizeof(PAIR_STATES[0]); ++i) {
                        printf("  %s\n", PAIR_STATES[i]);
                    }

                    exit(1);
                }

                break;
            }
            case KEYWORD_REDUCED_MASS:            input_block->reduced_mass = l->double_number; break;
            case KEYWORD_SO_POTENTIAL:            input_block->so_potential = strdup(l->string_storage.items); break;
            case KEYWORD_SO_DIPOLE:               input_block->so_dipole = strdup(l->string_storage.items); break;
            case KEYWORD_NITERATIONS:             params->niterations = l->int_number; break;
            case KEYWORD_TOTAL_TRAJECTORIES:      params->total_trajectories = l->int_number; break;
            case KEYWORD_CVODE_TOLERANCE:         params->cvode_tolerance = l->double_number; break;
            case KEYWORD_SAMPLING_TIME:           params->sampling_time = l->double_number; break;
            case KEYWORD_MAXTRAJECTORYLENGTH:     params->MaxTrajectoryLength = l->int_number; break;
            case KEYWORD_INITIALM0_NPOINTS:       params->initialM0_npoints = l->int_number; break;
            case KEYWORD_INITIALM2_NPOINTS:       params->initialM2_npoints = l->int_number; break;
            case KEYWORD_R0:                      params->R0 = l->double_number; break;
            case KEYWORD_APPROXIMATEFREQUENCYMAX: params->ApproximateFrequencyMax = l->double_number; break;
            case KEYWORD_TORQUE_CACHE_LEN:        params->torque_cache_len = l->int_number; break;
            case KEYWORD_TORQUE_LIMIT:            params->torque_limit = l->double_number; break;
            case KEYWORD_JINI_HISTOGRAM_BINS:     params->jini_histogram_bins = l->int_number; break;
            case KEYWORD_JINI_HISTOGRAM_MAX:      params->jini_histogram_max = l->double_number; break;
            case KEYWORD_JFIN_HISTOGRAM_BINS:     params->jfin_histogram_bins = l->int_number; break;
            case KEYWORD_JFIN_HISTOGRAM_MAX:      params->jfin_histogram_max = l->double_number; break;
            case KEYWORD_ORTHO_STATE_WEIGHT:      params->ortho_state_weight = l->double_number; break;
            case KEYWORD_PARA_STATE_WEIGHT:       params->para_state_weight = l->double_number; break;
            default: assert(false); 
        }
    }
}
            
void parse_monomer_block(Lexer *l, MonomerBlock *monomer_block) {
    while (true) {
        get_token(l);
        if ((l->token_type == TOKEN_BLOCK) && (strcasecmp(l->string_storage.items, "END") == 0)) {
            print_lexeme(l);
            return;
        }

        expect_token(l, TOKEN_KEYWORD);
        Keyword keyword_type = l->keyword_type;
        Token_Type expect_token = EXPECT_TOKEN[keyword_type]; 

        get_and_expect_token(l, TOKEN_EQ);
        get_and_expect_token(l, expect_token);

        switch (keyword_type) {
            case KEYWORD_MONOMER_TYPE: {
                if (strcasecmp(l->string_storage.items, "ATOM") == 0) {
                    monomer_block->t = ATOM;
                } else if (strcasecmp(l->string_storage.items, "LINEAR_MOLECULE") == 0) {
                    monomer_block->t = LINEAR_MOLECULE;
                } else if (strcasecmp(l->string_storage.items, "LINEAR_MOLECULE_REQ_INTEGER") == 0) {
                    monomer_block->t = LINEAR_MOLECULE_REQ_INTEGER;
                } else if (strcasecmp(l->string_storage.items, "LINEAR_MOLECULE_REQ_HALFINTEGER") == 0) {
                    monomer_block->t = LINEAR_MOLECULE_REQ_HALFINTEGER;
                } else if (strcasecmp(l->string_storage.items, "ROTOR") == 0) {
                    monomer_block->t = ROTOR;
                } else {
                    printf("ERROR: %s:%d:%d: unknown monomer type '%s'\n", l->loc.input_path, l->loc.line_number, l->loc.line_offset, l->string_storage.items);
                    printf("Available monomer_block types:\n");
                    for (size_t i = 0; i < sizeof(MONOMER_TYPES)/sizeof(MONOMER_TYPES[0]); ++i) {
                        printf(" %s\n", display_monomer_type(MONOMER_TYPES[i])); 
                    }
                    exit(1);
                }
                
                break;
            }
            case KEYWORD_DJ: monomer_block->DJ = l->double_number; break; 
            case KEYWORD_II: assert(false); break;
            default: assert(false);
        }

    }
} 

void parse_params(Lexer *l, CalcParams *params, InputBlock *input_block, MonomerBlock *monomer1_block, MonomerBlock *monomer2_block)
{
    while (true) { 
        get_token(l);
        if (is_eof(l)) return;

        expect_token(l, TOKEN_BLOCK);

        MonomerBlock *monomer_block = monomer1_block;

        if (strcasecmp(l->string_storage.items, "END") == 0) {
            printf("ERROR: %s:%d:%d: found '%s' without corresponding block beginning\n", l->loc.input_path, l->loc.line_number, l->loc.line_offset, l->string_storage.items);
            exit(1);
        } else if (strcasecmp(l->string_storage.items, "INPUT") == 0) {
            parse_input_block(l, input_block, params);
        } else if (strcasecmp(l->string_storage.items, "MONOMER") == 0) {
            parse_monomer_block(l, monomer_block); 
        }
    }
}

int main(int argc, char* argv[]) 
{
    UNUSED(argc);
    UNUSED(argv);

    const char *filename = "params.conf"; 
    String_Builder sb = {0};
    if (!read_entire_file(filename, &sb)) {
        exit(1);
    }

    Lexer l = lexer_new(filename, sb.items, sb.items + sb.count);
    //print_lexemes(&l);

    InputBlock input_block = {0}; 
    MonomerBlock monomer1  = {0};
    MonomerBlock monomer2  = {0};
    CalcParams params      = {0};

    parse_params(&l, &params,  &input_block, &monomer1, &monomer2);

    print_params(&params); 
    print_input_block(&input_block);
    print_monomer(&monomer1);

    return 0; 
}

/*
int main2(int argc, char* argv[])
{
    InputBlock input_block = {0};
    MonomerBlock monomer1  = {0};
    MonomerBlock monomer2  = {0};
    CalcParams params      = {0};
    parse_params_from_file("params.conf", &input_block, &monomer1, &monomer2, &params);
    
    print_input_block(&input_block);
    print_params(&params);


    // IDEA:
    // #embed directive в &INPUT, которая позволяет написать несколько строчек кода на С, который будет помещен во временный файл, скомпилирован в динамическую библиотеку
    // и мог быть вызван исполняемым файлом. 

    return 0;
}
*/
