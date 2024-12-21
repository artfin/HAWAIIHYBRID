#/bin/bash
set -xe

CC=gcc
CXX=g++
FLAGS="-Wall -Wextra -ggdb -O2 -march=native"

INC_SUNDIALS="-I/home/artfin/Desktop/lib/sundials-5.2.0/instdir/include"
INC_EIGEN="-I/usr/local/include/eigen3"

LIB_SUNDIALS="/home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_nvecserial.a /home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_cvode.a"

$CC $FLAGS -c -I./ ./PES-IDS/ai_pes_h2ar_leroy.c -o ./build/ai_pes_h2ar_leroy.o -lgsl -lgslcblas -lm 
$CC $FLAGS -c -I./ ./PES-IDS/ai_pes_co2ar.c -o ./build/ai_pes_co2ar.o -lgsl -lgslcblas -lm 
$CXX $FLAGS -c -I./ ./PES-IDS/ai_ids_co2ar.cpp -o ./build/ai_ids_co2ar.o -lgsl -lgslcblas -lm

$CC $FLAGS -c mtwist.c -o ./build/mtwist.o
$CC $FLAGS -c  -I./PES-IDS/ $INC_SUNDIALS hawaii.c -o ./build/hawaii.o
$CXX $FLAGS -c $INC_EIGEN angles_handler.cpp -o ./build/angles_handler.o

$CC $FLAGS -c array.c -o ./build/array.o
$CC $FLAGS $INC_SUNDIALS -c trajectory.c -o ./build/trajectory.o

SAMPLING_OBJ="./build/mtwist.o ./build/angles_handler.o ./build/hawaii.o ./build/array.o"
$CXX $FLAGS $INC_SUNDIALS $INC_EIGEN -I./ -I./PES-IDS/ $SAMPLING_OBJ ./build/ai_pes_co2ar.o ./build/ai_ids_co2ar.o ./examples/phase_space_integration_co2_ar.c -o ./examples/phase_space_integration_co2_ar.exe -lm $LIB_SUNDIALS -lgsl -lgslcblas -lstdc++  
$CXX $FLAGS $INC_SUNDIALS $INC_EIGEN -I./ -I./PES-IDS/ $SAMPLING_OBJ ./build/ai_pes_co2ar.o ./build/ai_ids_co2ar.o ./build/trajectory.o ./examples/trajectory_co2_ar.c -o ./examples/trajectory_co2_ar.exe -lm $LIB_SUNDIALS -lgsl -lgslcblas -lstdc++  
$CXX $FLAGS $INC_SUNDIALS $INC_EIGEN -I./ -I./PES-IDS/ $SAMPLING_OBJ ./build/trajectory.o ./build/ai_pes_h2ar_leroy.o ./examples/trajectory_h2_ar_requantized.c -o ./examples/trajectory_h2_ar_requantized.exe -lm $LIB_SUNDIALS -lgsl -lgslcblas -lstdc++  
