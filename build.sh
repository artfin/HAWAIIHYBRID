#/bin/bash
set -xe

FLAGS="-O2 -march=native"
INC="-I/home/artfin/Desktop/lib/sundials-5.2.0/instdir/include -I./PES-IDS/"
LIB="/home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_nvecserial.a /home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_cvode.a"
INC_EIGEN="-I/usr/local/include/eigen3"

g++ $FLAGS -c -I./ ./PES-IDS/ai_pes_co2ar.cpp -o ./build/ai_pes_co2ar.o -lgsl -lgslcblas -lm 

gcc $FLAGS -c mtwist.c -o ./build/mtwist.o
gcc $FLAGS -c $INC -Wall -Wextra hawaii.c -o ./build/hawaii.o
g++ $FLAGS -c $INC_EIGEN angles_handler.cpp -o ./build/angles_handler.o

g++ $FLAGS $INC $INC_EIGEN -ggdb -Wall -Wextra co2_ar.c ./build/ai_pes_co2ar.o ./build/hawaii.o ./build/mtwist.o ./build/angles_handler.o -o main.exe -lm $LIB -lgsl -lgslcblas -lstdc++  
