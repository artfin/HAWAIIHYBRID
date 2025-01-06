.PHONY: all clean test

CC    := gcc
CXX   := g++
FLAGS := -Wall -Wextra -ggdb -O2 -march=native

INC_SUNDIALS := -I/home/artfin/Desktop/lib/sundials-5.2.0/instdir/include
INC_EIGEN    := -I/usr/local/include/eigen3
LIB_GSL      := -lgsl -lgslcblas
LIB_SUNDIALS := /home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_nvecserial.a /home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_cvode.a

EXAMPLES := examples/phase_space_integration_co2_ar.exe \
			examples/trajectory_co2_ar.exe 				\
 			examples/trajectory_h2_ar_requantized.exe

all: $(EXAMPLES) 

clean: 
	rm -rf build/
	rm -v $(EXAMPLES) 

test: $(EXAMPLES)
	./run_tests.sh	


build/hawaii.o: hawaii.c hawaii.h | build
	$(CC) $(FLAGS) $(INC_SUNDIALS) -c $< -o $@ 

build/array.o: array.c array.h | build
	$(CC) $(FLAGS) -c $< -o $@ 

build/trajectory.o: trajectory.c trajectory.h | build
	$(CC) $(FLAGS) $(INC_SUNDIALS) -c $< -o $@ 

build/mtwist.o: mtwist.c mtwist.h | build
	$(CC) $(FLAGS) -c $< -o $@ 

build/ai_pes_co2_ar.o: ./PES-IDS/ai_pes_co2ar.c | build
	$(CC) $(FLAGS) -c -I./ $< -o $@ $(LINK_GSL) -lm 

build/ai_ids_co2_ar.o: ./PES-IDS/ai_ids_co2ar.cpp | build
	$(CC) $(FLAGS) -c -I./ $< -o $@ $(LINK_GSL) -lm 

build/ai_pes_h2ar_leroy.o: ./PES-IDS/ai_pes_h2ar_leroy.c | build
	$(CC) $(FLAGS) -c -I./ $< -o $@ $(LINK_GSL) -lm 

build/angles_handler.o: angles_handler.cpp angles_handler.hpp | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c $< -o $@ 

OBJ    := build/hawaii.o build/mtwist.o build/angles_handler.o build/array.o
CO2_AR := build/ai_pes_co2_ar.o build/ai_ids_co2_ar.o
H2_AR  := build/ai_pes_h2ar_leroy.o

examples/phase_space_integration_co2_ar.exe: examples/phase_space_integration_co2_ar.c build/hawaii.o $(OBJ) $(CO2_AR) 
	$(CXX) $(FLAGS) $(INC_SUNDIALS) $(INC_EIGEN) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/trajectory_co2_ar.exe: examples/trajectory_co2_ar.c build/trajectory.o $(OBJ) $(CO2_AR) 
	$(CXX) $(FLAGS) $(INC_SUNDIALS) $(INC_EIGEN) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/trajectory_h2_ar_requantized.exe: examples/trajectory_h2_ar_requantized.c build/trajectory.o $(OBJ) $(H2_AR)
	$(CXX) $(FLAGS) $(INC_SUNDIALS) $(INC_EIGEN) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  


build:
	mkdir -p $@

