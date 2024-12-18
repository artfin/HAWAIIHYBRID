#/bin/bash
set -xe

CC=gcc
CXX=g++
FLAGS="-O2 -march=native"

INC_SUNDIALS="-I/home/artfin/Desktop/lib/sundials-5.2.0/instdir/include"
INC_EIGEN="-I/usr/local/include/eigen3"

LIB_SUNDIALS="/home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_nvecserial.a /home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_cvode.a"

$CXX $FLAGS -c -I./ ./PES-IDS/ai_pes_co2ar.cpp -o ./build/ai_pes_co2ar.o -lgsl -lgslcblas -lm 
$CXX $FLAGS -c -I./ ./PES-IDS/ai_ids_co2ar.cpp -o ./build/ai_ids_co2ar.o -lgsl -lgslcblas -lm

$CC $FLAGS -c mtwist.c -o ./build/mtwist.o
$CC $FLAGS -c  -I./PES-IDS/ $INC_SUNDIALS -Wall -Wextra hawaii.c -o ./build/hawaii.o
$CXX $FLAGS -c $INC_EIGEN angles_handler.cpp -o ./build/angles_handler.o

$CC $FLAGS -c array.c -o ./build/array.o
$CC $FLAGS $INC_SUNDIALS -c trajectory.c -o ./build/trajectory.o

SAMPLING_OBJ="./build/ai_pes_co2ar.o ./build/ai_ids_co2ar.o ./build/mtwist.o ./build/angles_handler.o ./build/hawaii.o ./build/array.o"
$CXX $FLAGS $INC_SUNDIALS $INC_EIGEN -I./ -I./PES-IDS/ -ggdb -Wall -Wextra $SAMPLING_OBJ ./examples/phase_space_integration_co2_ar.c -o ./examples/phase_space_integration_co2_ar.exe -lm $LIB_SUNDIALS -lgsl -lgslcblas -lstdc++  
$CXX $FLAGS $INC_SUNDIALS $INC_EIGEN -I./ -I./PES-IDS/ -ggdb -Wall -Wextra ./examples/trajectory_co2_ar.c ./build/ai_pes_co2ar.o ./build/hawaii.o ./build/mtwist.o ./build/angles_handler.o ./build/trajectory.o ./build/array.o -o ./examples/trajectory_co2_ar.exe -lm $LIB_SUNDIALS -lgsl -lgslcblas -lstdc++  
