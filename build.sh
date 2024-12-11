#/bin/bash
set -xe

INC="-I/home/artfin/Desktop/lib/sundials-5.2.0/instdir/include"
LIB="/home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_nvecserial.a /home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_cvode.a"

gcc $INC -ggdb -Wall -Wextra main.c -o main.exe -lm $LIB 
