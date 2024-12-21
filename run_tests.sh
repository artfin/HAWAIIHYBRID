#!/bin/bash

set -xe

TESTS=("phase_space_integration_co2_ar.exe" "trajectory_co2_ar.exe")

for test in "${TESTS[@]}"; do
    output_file="${test%.exe}.out"
    if ! "./examples/$test" > "./examples/$output_file" 2>&1; then
        echo "ERROR: $test returned non-zero exit code"
        exit 1
    fi
done

echo "All tests ran successfully"
