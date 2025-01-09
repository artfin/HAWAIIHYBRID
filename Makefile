.PHONY: all clean test docs

CC     := gcc
CXX    := g++
MPICC  := mpicc
MPICXX := mpic++
FLAGS  := -Wall -Wextra -ggdb -O2 -march=native

INC_SUNDIALS := -I/home/artfin/Desktop/lib/sundials-5.2.0/instdir/include
INC_EIGEN    := -I/usr/local/include/eigen3
INC_HEP      := -I/home/artfin/Desktop/lib/hep-mc-0.7/include/
INC          := $(INC_SUNDIALS) $(INC_EIGEN) $(INC_HEP) 
LIB_GSL      := -lgsl -lgslcblas
LIB_SUNDIALS := /home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_nvecserial.a /home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_cvode.a

EXAMPLES := examples/phase_space_integration_co2_ar.exe \
			examples/trajectory_co2_ar.exe 				\
 			examples/trajectory_h2_ar_requantized.exe   \
			examples/correlation_co2_ar.exe

all: $(EXAMPLES) 

clean: 
	rm -rf -v build/
	rm -v $(EXAMPLES)
	rm -v docs/main.toc
	rm -v docs/main.aux
	rm -v docs/main.log

test: $(EXAMPLES)
	./run_tests.sh	

docs: docs/main.tex
	cd docs && pdflatex main.tex && pdflatex main.tex

build/hawaii.o: hawaii.c | build
	$(CC) $(FLAGS) $(INC) -c -MD $< -o $@ 

build/mpi_hawaii.o: hawaii.c | build
	$(MPICC) $(FLAGS) $(INC) -DUSE_MPI -c -MD $< -o $@ 

build/hep_hawaii.o: hep_hawaii.cpp | build
	$(MPICC) $(FLAGS) $(INC) -c -MD $< -o $@ 

build/array.o: array.c | build
	$(CC) $(FLAGS) -c -MD $< -o $@ 

build/trajectory.o: trajectory.c | build
	$(CC) $(FLAGS) $(INC) -c -MD $< -o $@ 

build/mtwist.o: mtwist.c | build
	$(CC) $(FLAGS) -c -MD $< -o $@ 

build/ai_pes_co2_ar.o: ./PES-IDS/ai_pes_co2ar.c | build
	$(CC) $(FLAGS) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 

build/ai_ids_co2_ar.o: ./PES-IDS/ai_ids_co2ar.cpp | build
	$(CC) $(FLAGS) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 

build/ai_pes_h2ar_leroy.o: ./PES-IDS/ai_pes_h2ar_leroy.c | build
	$(CC) $(FLAGS) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 

build/angles_handler.o: angles_handler.cpp | build
	$(CXX) $(FLAGS) $(INC) -c -MD $< -o $@ 

OBJ     := build/hawaii.o build/mtwist.o build/angles_handler.o build/array.o build/trajectory.o
MPI_OBJ := build/mpi_hawaii.o build/mtwist.o build/angles_handler.o build/array.o build/trajectory.o
CO2_AR  := build/ai_pes_co2_ar.o build/ai_ids_co2_ar.o
H2_AR   := build/ai_pes_h2ar_leroy.o

examples/phase_space_integration_co2_ar.exe: examples/phase_space_integration_co2_ar.cpp build/hawaii.o $(OBJ) $(CO2_AR) 
	$(CXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/mpi_phase_space_integration_co2_ar.exe: examples/mpi_phase_space_integration_co2_ar.cpp $(MPI_OBJ) build/hep_hawaii.o $(CO2_AR)  
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/trajectory_co2_ar.exe: examples/trajectory_co2_ar.c build/trajectory.o $(OBJ) $(CO2_AR) 
	$(CXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/trajectory_h2_ar_requantized.exe: examples/trajectory_h2_ar_requantized.c build/trajectory.o $(OBJ) $(H2_AR)
	$(CXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/correlation_co2_ar.exe: examples/correlation_co2_ar.cpp build/trajectory.o $(MPI_OBJ) $(CO2_AR) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) 


build:
	mkdir -p $@

include $(wildcard build/*.d)
