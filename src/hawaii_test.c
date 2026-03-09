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


Report EXPECTED_TESTS_STATUS[] = {
    { .name = "empty-processing.conf",                         .run_status = Success },
    { .name = "read-cf.conf",                                  .run_status = Fail },
    { .name = "read-sf.conf",                                  .run_status = Fail },
    { .name = "convert-cf-to-sf.conf",                         .run_status = Success },
    { .name = "cf-to-sf-argtype.conf",                         .run_status = Fail },
    { .name = "string-literal.conf",                           .run_status = Fail },
    { .name = "naked-drop.conf",                               .run_status = Fail },
    { .name = "naked-dup.conf",                                .run_status = Fail },
    { .name = "cmp.conf",                                      .run_status = Fail },
    { .name = "fit-baseline-nargs.conf",                       .run_status = Fail },
    { .name = "average-cfs.conf",                              .run_status = Fail },
    { .name = "compute-classical-moments.conf",                .run_status = Success },
    { .name = "compute-classical-moments-with-truncation.conf",.run_status = Success },
    { .name = "compute-classical-moments-nargs.conf",          .run_status = Fail },
    { .name = "compute-classical-moments-argname.conf",        .run_status = Fail },
    { .name = "compute-classical-moments-argtype.conf",        .run_status = Fail },
    { .name = "compute-quantum-moments.conf",                  .run_status = Fail },
    { .name = "compute-quantum-moments-with-truncation.conf",  .run_status = Success },
    { .name = "compute-quantum-moments-nargs.conf",            .run_status = Fail },
    { .name = "compute-quantum-moments-argname.conf",          .run_status = Fail },
    { .name = "compute-quantum-moments-argtype.conf",          .run_status = Fail },
    { .name = "compute-alpha.conf",                            .run_status = Fail },
    { .name = "desymmetrizations.conf",                        .run_status = Success },
    { .name = "d1-argtype-sp.conf",                            .run_status = Fail },
    { .name = "write-cf.conf",                                 .run_status = Success },
    { .name = "write-sf.conf",                                 .run_status = Success },
    { .name = "write-sp.conf",                                 .run_status = Success },
    { .name = "add-spectra.conf",                              .run_status = Success },
    { .name = "dup2.conf",                                     .run_status = Fail },
    { .name = "drop2.conf",                                    .run_status = Fail },
    { .name = "smooth.conf",                                   .run_status = Success },
    { .name = "inv-d3.conf",                                   .run_status = Success },
    { .name = "inv-d2.conf",                                   .run_status = Success },
    { .name = "swap.conf",                                     .run_status = Success },
    { .name = "rot.conf",                                      .run_status = Success },
};

#define TEST_COUNT sizeof(EXPECTED_TESTS_STATUS)/sizeof(EXPECTED_TESTS_STATUS[0]) 
static_assert(TEST_COUNT == 35, "");

Status run_test(Cmd *cmd, const char *test_name) {
    cmd_append(cmd, "./driver.exe", "-quiet", temp_sprintf("./tests/%s", test_name));
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
        if (strcmp(EXPECTED_TESTS_STATUS[j].name, test_name) == 0) return true;
    }

    printf("ERROR: test '%s' is not found\n", test_name);
    return false; 
}


Status expected_run_status(const char *test_name) {
    for (size_t j = 0; j < TEST_COUNT; ++j) {
        if (strcmp(EXPECTED_TESTS_STATUS[j].name, test_name) == 0)
            return EXPECTED_TESTS_STATUS[j].run_status;
    }
    return Success;
}

void print_reports(Reports *reports) {
    printf("Reports:\n");
    printf("\t%-40s\t%-10s\t%s\n", "Test Name", "Run", "Result");
    for (size_t i = 0; i < reports->count; ++i) {
        Report *report = &reports->items[i];
        Status expected_run = expected_run_status(report->name);

        printf("%2zu \t%-40s\t", i+1, report->name);

        if (report->run_status == expected_run) {
            printf("\e[32m%-10s\e[0m\t", STATUS_AS_STR[report->run_status]);
        } else {
            printf("\e[31m%-10s\e[0m\t", STATUS_AS_STR[report->run_status]);
        }

        if (report->result_status == Success) {
            printf("\e[32m%s\e[0m", STATUS_AS_STR[report->result_status]);
        } else {
            printf("\e[31m%s\e[0m", STATUS_AS_STR[report->result_status]);
        }

        printf("\n");
    }
}

void collect_test_reports(Reports *reports)
{
    Cmd cmd = {0};
    String_Builder tmp_filename_content = {0};
    String_Builder out_filename_content = {0};

    for (size_t i = 0; i < TEST_COUNT; ++i) {
        const char *test_name = EXPECTED_TESTS_STATUS[i].name;
        Report report = {
            .name = strdup(test_name),
            .run_status = run_test(&cmd, test_name),
        };

        const char *tmp_filename = temp_sprintf(TMP_FILENAME_TEMPLATE, test_name);
        if (!read_entire_file(tmp_filename, &tmp_filename_content)) {
            printf("ERROR: could not read the file %s\n", tmp_filename);
            continue;
        }
        sb_append_null(&tmp_filename_content);

        const char *out_filename = temp_sprintf(OUT_FILENAME_TEMPLATE, test_name);
        if (file_exists(out_filename)) {
            if (!read_entire_file(out_filename, &out_filename_content)) {
                printf("ERROR: could not read the file %s\n", out_filename);
                continue;
            }
            sb_append_null(&out_filename_content);
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

void cmd_replay_all(void) {
    Reports reports = {0};
    collect_test_reports(&reports);
    print_reports(&reports);
}

void cmd_replay_selected(int rest_argc, char **rest_argv) {
    Cmd cmd = {0};
    String_Builder tmp_filename_content = {0};
    String_Builder out_filename_content = {0};
    Reports reports = {0};

    for (int i = 0; i < rest_argc; ++i) {
        const char *test_name = *rest_argv;
        rest_argv += 1;

        if (!test_exists(test_name)) continue;

        Status run_status = run_test(&cmd, test_name);

        const char *tmp_filename = temp_sprintf(TMP_FILENAME_TEMPLATE, test_name);
        if (!read_entire_file(tmp_filename, &tmp_filename_content)) {
            printf("ERROR: could not read the file %s\n", tmp_filename);
            continue;
        }
        sb_append_null(&tmp_filename_content);

        printf("%s", tmp_filename_content.items);

        Report report = {
            .name = strdup(test_name),
            .run_status = run_status,
        };

        const char *out_filename = temp_sprintf(OUT_FILENAME_TEMPLATE, test_name);
        if (file_exists(out_filename)) {
            if (!read_entire_file(out_filename, &out_filename_content)) {
                printf("ERROR: could not read the file %s\n", out_filename);
                continue;
            }
            sb_append_null(&out_filename_content);
            report.result_status = (strcmp(tmp_filename_content.items, out_filename_content.items) == 0) ? Success : Fail;
            out_filename_content.count = 0;
        } else {
            printf("ERROR: missing expected output file for '%s'\n", test_name);
        }

        da_append(&reports, report);
        tmp_filename_content.count = 0;
    }

    sb_free(tmp_filename_content);
    sb_free(out_filename_content);

    print_reports(&reports);
}

void cmd_record(int rest_argc, char **rest_argv) {
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

void cmd_diff(int rest_argc, char **rest_argv) {
    if (rest_argc == 0) {
        printf("ERROR: no test names are provided to diff\n");
        exit(1);
    }

    Cmd cmd = {0};

    for (int i = 0; i < rest_argc; ++i) {
        const char *test_name = *rest_argv;
        rest_argv += 1;

        if (!test_exists(test_name)) continue;

        run_test(&cmd, test_name);

        const char *out_filename = temp_sprintf(OUT_FILENAME_TEMPLATE, test_name);
        const char *tmp_filename = temp_sprintf(TMP_FILENAME_TEMPLATE, test_name);

        if (!file_exists(out_filename)) {
            printf("ERROR: missing recorded output for '%s'\n", test_name);
            continue;
        }

        cmd_append(&cmd, "diff", "--color", out_filename, tmp_filename);
        nob_cmd_run_sync_and_reset(&cmd);
    }
}

int main(int argc, char *argv[])
{
    bool *replay = flag_bool("replay", false, "Replay all tests (if no argument is provided) or a particular test");
    bool *record = flag_bool("record", false, "Record the result for particular test");
    bool *diff = flag_bool("diff", false, "Run a test and show diff between recorded and obtained output");
    bool *help = flag_bool("help", false, "Print this help message");

    if (!flag_parse(argc, argv)) {
        usage();
        exit(1);
    }

    if (*help) {
        usage();
        return 0;
    }

    int rest_argc = flag_rest_argc();
    char **rest_argv = flag_rest_argv();

    if (*replay) {
        if (rest_argc == 0) {
            cmd_replay_all();
        } else {
            cmd_replay_selected(rest_argc, rest_argv);
        }
    }

    if (*record) {
        cmd_record(rest_argc, rest_argv);
    }

    if (*diff) {
        cmd_diff(rest_argc, rest_argv);
    }

    return 0;
}

/*
 *  Copyright (C) 2026 A.Finenko & D.Chistikov
 *  Distributed under the GNU General Public License, version 3
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */       

