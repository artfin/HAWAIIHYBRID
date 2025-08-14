#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>

#define NOB_STRIP_PREFIX
#define NOB_IMPLEMENTATION
#include "nob.h"

#define FLAG_IMPLEMENTATION
#include "flag.h"

typedef enum {
    Success,
    Fail,
} Status;

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
    "parse-cf-header.conf",
};
#define TEST_COUNT sizeof(TEST_NAMES)/sizeof(TEST_NAMES[0]) 

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

int main(int argc, char *argv[])
{
    bool *replay = flag_bool("replay", false, "Replay all tests or a particular test");
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

        Cmd cmd = {0};

        if (rest_argc == 0) {
            Reports reports = {0};
                
            String_Builder tmp_filename_content = {0};
            String_Builder out_filename_content = {0};

            for (size_t i = 0; i < TEST_COUNT; ++i) {
                const char *test_name = TEST_NAMES[i];
                Report report = {
                    .name = strdup(test_name),
                    .run_status = run_test(&cmd, test_name),
                };

                const char *tmp_filename = temp_sprintf("./tests/%s.out.tmp", test_name);
                if (!read_entire_file(tmp_filename, &tmp_filename_content)) {
                    printf("ERROR: could not read the file %s\n", tmp_filename);
                    continue;     
                }
                sb_append_null(&tmp_filename_content); // just in case 
                 
                const char *out_filename = temp_sprintf("./tests/%s.out", test_name);
                if (file_exists(out_filename)) {
                    if (!read_entire_file(out_filename, &out_filename_content)) {
                        printf("ERROR: could not read the file %s\n", out_filename);
                        continue;
                    }
                }
                sb_append_null(&out_filename_content); // just in case 
                
                report.result_status = (strcmp(tmp_filename_content.items, out_filename_content.items) == 0) ? Success : Fail; 
                da_append(&reports, report);
                
                tmp_filename_content.count = 0;
                out_filename_content.count = 0;
            }

            sb_free(tmp_filename_content);
            sb_free(out_filename_content);

            printf("Reports:\n");
            printf("\t test name \t\t run \t result\n"); 
            for (size_t i = 0; i < reports.count; ++i) {
                Report *report = &reports.items[i];
                printf("\t%s \t", report->name);
                if (report->run_status == Success) {
                    printf("\e[32mSuccess\e[0m");
                } else if (report->run_status == Fail) {
                    printf("\e[31mFail\e[0m");
                }

                if (report->result_status == Success) {
                    printf("\t \e[32mSuccess\e[0m");
                } else if (report->result_status == Fail) {
                    printf("\t \e[31mFail\e[0m");
                }

                printf("\n");
            }
        } else {
            for (int i = 0; i < rest_argc; ++i) {
                const char *test_name = *rest_argv;
                rest_argv += 1;

                bool found = false; 
                for (size_t j = 0; j < TEST_COUNT; ++j) {
                    if (strcmp(TEST_NAMES[j], test_name) == 0) {
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    printf("ERROR: test %s is not found\n", test_name);
                    return 1; 
                }
                
                Status status = run_test(&cmd, test_name);
                (void) status;
            }
        }
    }

    if (*record) {
        int rest_argc = flag_rest_argc();
        char **rest_argv = flag_rest_argv();

        (void) rest_argc;
        (void) rest_argv; 
    }

    return 0;
}
