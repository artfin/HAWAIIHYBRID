#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>

#define NOB_STRIP_PREFIX
#define NOB_IMPLEMENTATION
#include "nob.h"

#define FLAG_IMPLEMENTATION
#include "flag.h"

// the result of test recording will be saved using this template
#define OUT_FILENAME_TEMPLATE "./tests/%s.out"
// the temporary result when a test is just replayed is saved using this template
#define TMP_FILENAME_TEMPLATE "./tests/%s.out.tmp"

// TODO: allow stripping off '.conf' extension when specifying the test name

typedef enum {
    Success,
    Fail,
} Status;

const char *STATUS_AS_STR[] = {
    [Success] = "Success",
    [Fail]    = "Fail",
};

typedef struct {
    const char *name;
    Status run_status; // whether running was successful
    Status result_status; // whether the result of running was correct 
} Report;

typedef struct {
    Report *items;
    size_t count;
    size_t capacity;
} Reports;

const char *TEST_NAMES[] = {
    "empty-processing.conf",
    "read-cf.conf",
    "read-sf.conf",
    "convert-cf-to-sf.conf",
    "string-literal.conf",
    "naked-drop.conf",
    "naked-dup.conf",
    "cmp.conf",
    "fit-baseline-nargs.conf",
    "average-cfs.conf",
    "compute-classical-moments.conf",
    "compute-classical-moments-with-truncation.conf",
    "compute-classical-moments-nargs.conf",
    "compute-classical-moments-argname.conf",
    "compute-classical-moments-argtype.conf",
    "compute-quantum-moments.conf",
    "compute-quantum-moments-with-truncation.conf",
    "compute-quantum-moments-nargs.conf",
    "compute-quantum-moments-argname.conf",
    "compute-quantum-moments-argtype.conf",
    "compute-alpha.conf",
};
#define TEST_COUNT sizeof(TEST_NAMES)/sizeof(TEST_NAMES[0]) 

Status EXPECTED_RUN_STATUS[TEST_COUNT] = {
    Success,
    Fail,
    Fail,
    Success,
    Fail,
    Fail,
    Fail,
    Fail,
    Fail,
    Fail,
    Success,
    Success,
    Fail,
    Fail,
    Fail,
    Success,
    Success,
    Fail,
    Fail,
    Fail,
    Fail,
};

static_assert(TEST_COUNT == 21, "");

Status run_test(Cmd *cmd, const char *test_name) {
    cmd_append(cmd, "./driver.exe", temp_sprintf("./tests/%s", test_name));
    Fd fdout = fd_open_for_write(temp_sprintf("./tests/%s.out.tmp", test_name));
    if (!nob_cmd_run_sync_redirect_and_reset(cmd, (Cmd_Redirect) { .fdout = &fdout })) {
        return Fail;
    }

    return Success;
}

void usage(void)
{
    printf("Hawaii Hybrid Testing Tool\n");
    printf("Usage: %s [OPTIONS]\n", flag_program_name());
    printf("  OPTIONS:\n");
    flag_print_options(stdout); 
}


bool test_exists(const char *test_name) 
{
    for (size_t j = 0; j < TEST_COUNT; ++j) {
        if (strcmp(TEST_NAMES[j], test_name) == 0) return true;
    }

    printf("ERROR: test '%s' is not found\n", test_name);
    return false; 
}


void collect_test_reports(Reports *reports) 
{
    Cmd cmd = {0};
    String_Builder tmp_filename_content = {0};
    String_Builder out_filename_content = {0};

    for (size_t i = 0; i < TEST_COUNT; ++i) {
        const char *test_name = TEST_NAMES[i];
        Report report = {
            .name = strdup(test_name),
            .run_status = run_test(&cmd, test_name),
        };

        const char *tmp_filename = temp_sprintf(TMP_FILENAME_TEMPLATE, test_name);
        if (!read_entire_file(tmp_filename, &tmp_filename_content)) {
            printf("ERROR: could not read the file %s\n", tmp_filename);
            continue;     
        }
        sb_append_null(&tmp_filename_content); // just in case 

        const char *out_filename = temp_sprintf(OUT_FILENAME_TEMPLATE, test_name);
        if (file_exists(out_filename)) {
            if (!read_entire_file(out_filename, &out_filename_content)) {
                printf("ERROR: could not read the file %s\n", out_filename);
                continue;
            }

            sb_append_null(&out_filename_content); // just in case 
            
            report.result_status = (strcmp(tmp_filename_content.items, out_filename_content.items) == 0) ? Success : Fail; 
            da_append(reports, report);
        
            out_filename_content.count = 0;
        } else {
            printf("ERROR: missing expected output file for '%s'\n", test_name);
        }

        tmp_filename_content.count = 0;
    }

    sb_free(tmp_filename_content);
    sb_free(out_filename_content);
}

int main(int argc, char *argv[])
{
    bool *replay = flag_bool("replay", false, "Replay all tests (if no argument is provided) or a particular test");
    bool *record = flag_bool("record", false, "Record the result for particular test");
    bool *help = flag_bool("help", false, "Print this help message");

    if (!flag_parse(argc, argv)) {
        usage();
        exit(1);
    }

    if (*help) {
        usage();
        return 0;
    }

    if (*replay) {
        int rest_argc = flag_rest_argc();
        char **rest_argv = flag_rest_argv();


        if (rest_argc == 0) {
            Reports reports = {0};
            collect_test_reports(&reports);

            printf("Reports:\n");
            printf("\t%-40s\t%-10s\t%s\n", "Test Name", "Run", "Result"); 
            for (size_t i = 0; i < reports.count; ++i) {
                Report *report = &reports.items[i];
                printf("%2zu \t%-40s\t", i+1, report->name);

                if (report->run_status == EXPECTED_RUN_STATUS[i]) {
                    printf("\e[32m%-10s\e[0m\t", STATUS_AS_STR[report->run_status]);
                } else {
                    printf("\e[31m%-10s\e[0m\t", STATUS_AS_STR[report->run_status]);
                }

                if (report->result_status == Success) {
                    printf("\e[32m%s\e[0m", STATUS_AS_STR[report->result_status]);
                } else if (report->result_status == Fail) {
                    printf("\e[31m%s\e[0m", STATUS_AS_STR[report->result_status]);
                }

                printf("\n");
            }
        } else {
            Cmd cmd = {0};

            String_Builder tmp_filename_content = {0};

            for (int i = 0; i < rest_argc; ++i) {
                const char *test_name = *rest_argv;
                rest_argv += 1;

                if (!test_exists(test_name)) continue; 
                
                Status status = run_test(&cmd, test_name);
                (void) status;

                const char *tmp_filename = temp_sprintf(TMP_FILENAME_TEMPLATE, test_name);
                if (!read_entire_file(tmp_filename, &tmp_filename_content)) {
                    printf("ERROR: could not read the file %s\n", tmp_filename);
                    continue;     
                }
                sb_append_null(&tmp_filename_content); // just in case 

                printf("%s\n\n", tmp_filename_content.items);
            }

            sb_free(tmp_filename_content);
        }
    }

    if (*record) {
        int rest_argc = flag_rest_argc();
        char **rest_argv = flag_rest_argv();

        if (rest_argc == 0) {
            printf("ERROR: no test names are provided to record\n");
            exit(1);
        }

        Cmd cmd = {0};
        
        for (int i = 0; i < rest_argc; ++i) {
            const char *test_name = *rest_argv;
            rest_argv += 1;

            if (!test_exists(test_name)) continue;

            printf("Recording test '%s'\n", test_name);
            Status status = run_test(&cmd, test_name);
            (void) status;

            cmd_append(&cmd, "mv", temp_sprintf("./tests/%s.out.tmp", test_name), temp_sprintf("./tests/%s.out", test_name));
            if (!nob_cmd_run_sync_and_reset(&cmd)) {
                printf("ERROR: failed to record test '%s'\n", test_name);
                continue;
            }
        }
    }

    return 0;
}
