.PHONY: all clean test docs

CC     := gcc
F      := gfortran
CXX    := g++
MPICC  := mpicc
MPICXX := mpic++

# -Wswitch-enum: if default statement is present in the switch case, but not all the enum values are covered, the warning will still be emitted 
FLAGS_DEBUG   := -Wall -Wextra -Wswitch-enum -ggdb -O0 
FLAGS_RELEASE := -Wall -Wextra -Wswitch-enum -O2 -march=native -mtune=native # -pg -ggdb
FLAGS_EIGEN   := -Wall -Wextra -Wswitch-enum -O2 # -pg -ggdb 
FLAGS := $(FLAGS_DEBUG)

INC_SUNDIALS := -I/home/artfin/Desktop/lib/sundials-5.2.0/instdir/include
INC_EIGEN    := -I/usr/local/include/eigen3
INC_HEP      := -I/home/artfin/Desktop/lib/hep-mc-0.7/include/
INC          := $(INC_SUNDIALS) $(INC_EIGEN) $(INC_HEP) 
LIB_GSL      := -lgsl -lgslcblas
LIB_SUNDIALS := /home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_nvecserial.a /home/artfin/Desktop/lib/sundials-5.2.0/instdir/lib/libsundials_cvode.a

EXAMPLES := examples/phase_space_integration_co2_ar.exe      \
			examples/mpi_phase_space_integration_co2_ar.exe  \
			examples/mpi_phase_space_integration_ch4_co2.exe \
			examples/phase_space_integration_he_ar.exe       \
			examples/trajectory_co2_ar.exe 				     \
 			examples/trajectory_h2_ar_requantized.exe        \
			examples/trajectory_ch4_co2.exe 		         \
		    examples/correlation_he_ar.exe                   \
			examples/correlation_co_ar.exe                   \
			examples/correlation_co2_ar.exe                  \
			examples/correlation_n2_ar.exe                   \
			examples/correlation_array_n2_ar.exe             \
			examples/correlation_array_co2_ar.exe            \
			examples/correlation_ch4_co2.exe                 \
			examples/correlation_array_ch4_co2.exe           \
			examples/prmu_calculation_co2_ar.exe             \
			examples/prmu_calculation_co_ar_requantized.exe  \
			examples/prmu_calculation_h2_ar_requantized.exe  \
			examples/fftrump.exe                             \
			examples/test_sb.exe                             \
			examples/test_loess.exe                          \
			examples/test_fft.exe

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
	cd docs && pdflatex main.tex && bibtex main && pdflatex main.tex

build/hawaii.o: hawaii.c | build
	$(CC) $(FLAGS) $(INC) -c -MD $< -o $@ 

build/mpi_hawaii.o: hawaii.c | build
	$(MPICC) $(FLAGS) $(INC) -DUSE_MPI -c -MD $< -o $@ 

build/hep_hawaii.o: hep_hawaii.cpp | build
	$(MPICXX) $(FLAGS) $(INC) -c -MD $< -o $@ 

build/array.o: array.c | build
	$(CC) $(FLAGS) -c -MD $< -o $@ 

build/trajectory.o: trajectory.c | build
	$(CC) $(FLAGS) $(INC) -c -MD $< -o $@ 

build/mtwist.o: mtwist.c | build
	$(CC) $(FLAGS) -c -MD $< -o $@ 

build/angles_handler.o: angles_handler.cpp | build
	$(CXX) $(FLAGS_EIGEN) $(INC) -c -MD $< -o $@ 

build/loess.o: loess.cpp | build
	$(CXX) $(FLAGS_EIGEN) $(INC) -c -MD $< -o $@ -fopenmp  

###########################################################
###################### CO2-Ar #############################
###########################################################
build/ai_pes_co2_ar.o: ./PES-IDS/ai_pes_co2ar.c | build
	$(CC) $(FLAGS) $(INC) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 

build/ai_ids_co2_ar.o: ./PES-IDS/ai_ids_co2ar.cpp | build
	$(CC) $(FLAGS) $(INC) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 
###########################################################

###########################################################
###################### H2-Ar ##############################
###########################################################
build/c_basis_2_2_1_3_intermolecular.o: ./PES-IDS/c_basis_2_2_1_3_intermolecular.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -I./ $< -o $@ -lm

build/c_basis_2_1_1_1_3_intermolecular.o: ./PES-IDS/c_basis_2_1_1_1_3_intermolecular.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -I./ $< -o $@ -lm

build/c_basis_1_1_2_1_3_intermolecular.o: ./PES-IDS/c_basis_1_1_2_1_3_intermolecular.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -I./ $< -o $@ -lm

build/ai_pes_h2ar_leroy.o: ./PES-IDS/ai_pes_h2ar_leroy.c | build
	$(CC) $(FLAGS) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm

build/ai_ids_h2_ar_pip_nn.o: ./PES-IDS/ai_ids_h2_ar_pip_nn.cpp | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c -MD -I./ $< -o $@ -lm
###########################################################

###########################################################
###################### CO-Ar ##############################
###########################################################
build/potv.o: ./PES-IDS/potv.f | build
	$(F) -c $< -o $@ 

build/potv_d.o: ./PES-IDS/potv_d.f03 | build
	$(F) -c $< -o $@
###########################################################

###########################################################
###################### N2-Ar ##############################
###########################################################
build/c_basis_2_1_4_purify.o: ./PES-IDS/c_basis_2_1_4_purify.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -I./ $< -o $@ -lm
   
build/c_jac_2_1_4_purify.o: ./PES-IDS/c_jac_2_1_4_purify.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -I./ $< -o $@ -lm

build/c_basis_2_2_1_3_purify.o: ./PES-IDS/c_basis_2_2_1_3_purify.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -I./ $< -o $@ -lm

build/c_basis_2_1_1_1_3_purify.o: ./PES-IDS/c_basis_2_1_1_1_3_purify.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -I./ $< -o $@ -lm

build/c_basis_1_1_2_1_3_purify.o: ./PES-IDS/c_basis_1_1_2_1_3_purify.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -I./ $< -o $@ -lm

build/cnpy.o: ./PES-IDS/cnpy.cpp | build
	$(CC) $(FLAGS) -c -MD ./ $< -o $@ -lm

build/ai_pes_n2_ar_pip_nn.o: ./PES-IDS/ai_pes_n2_ar_pip_nn.cpp | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c -MD -I./ $< -o $@ -lm

build/ai_ids_n2_ar_pip_nn.o: ./PES-IDS/ai_ids_n2_ar_pip_nn.cpp | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c -MD -I./ $< -o $@ -lm
###########################################################

###########################################################
##################### CH4-CO2 #############################
###########################################################
# build/ai_pes_ch4_co2.o: ./PES-IDS/ai_pes_ch4_co2.c | build
# 	$(CC) $(FLAGS) $(INC_GSL) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 
# 
# build/ai_pes_ch4_co2_dEdR.o: ./PES-IDS/ai_pes_ch4_co2_dEdR.c | build
# 	$(CC) $(FLAGS) $(INC_GSL) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 
# 
# build/ai_pes_ch4_co2_dEdphi1.o: ./PES-IDS/ai_pes_ch4_co2_dEdphi1.c | build
# 	$(CC) $(FLAGS) $(INC_GSL) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 
# 
# build/ai_pes_ch4_co2_dEdtheta1.o: ./PES-IDS/ai_pes_ch4_co2_dEdtheta1.c | build
# 	$(CC) $(FLAGS) $(INC_GSL) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 
# 
# build/ai_pes_ch4_co2_dEdphi2.o: ./PES-IDS/ai_pes_ch4_co2_dEdphi2.c | build
# 	$(CC) $(FLAGS) $(INC_GSL) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 
# 
# build/ai_pes_ch4_co2_dEdtheta2.o: ./PES-IDS/ai_pes_ch4_co2_dEdtheta2.c | build
# 	$(CC) $(FLAGS) $(INC_GSL) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 
# 
# build/ai_ids_ch4_co2.o: ./PES-IDS/ai_ids_ch4_co2.cpp | build
# 	$(CC) $(FLAGS) $(INC_GSL) $(INC_EIGEN) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 
###########################################################

OBJ     := build/hawaii.o build/mtwist.o build/angles_handler.o build/array.o build/trajectory.o
MPI_OBJ := build/mpi_hawaii.o build/mtwist.o build/angles_handler.o build/array.o build/trajectory.o
CO2_AR  := build/ai_pes_co2_ar.o build/ai_ids_co2_ar.o
N2_AR   := build/cnpy.o -lz build/ai_pes_n2_ar_pip_nn.o build/ai_ids_n2_ar_pip_nn.o \
		   build/c_basis_2_1_4_purify.o build/c_jac_2_1_4_purify.o \
		   build/c_basis_2_2_1_3_purify.o build/c_basis_2_1_1_1_3_purify.o build/c_basis_1_1_2_1_3_purify.o
H2_AR   := build/cnpy.o -lz build/ai_pes_h2ar_leroy.o build/ai_ids_h2_ar_pip_nn.o \
		   build/c_basis_2_2_1_3_intermolecular.o build/c_basis_2_1_1_1_3_intermolecular.o build/c_basis_1_1_2_1_3_intermolecular.o
CO_AR   := build/potv.o build/potv_d.o 
CH4_CO2 := build/ai_pes_ch4_co2.o build/ai_pes_ch4_co2_dEdR.o build/ai_pes_ch4_co2_dEdphi1.o build/ai_pes_ch4_co2_dEdtheta1.o \
		   build/ai_pes_ch4_co2_dEdphi2.o build/ai_pes_ch4_co2_dEdtheta2.o build/ai_ids_ch4_co2.o
 
 

examples/phase_space_integration_co2_ar.exe: examples/phase_space_integration_co2_ar.cpp build/hawaii.o $(OBJ) $(CO2_AR) 
	$(CXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/phase_space_integration_he_ar.exe: examples/phase_space_integration_he_ar.cpp build/hep_hawaii.o $(MPI_OBJ) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/mpi_phase_space_integration_co2_ar.exe: examples/mpi_phase_space_integration_co2_ar.cpp $(MPI_OBJ) build/hep_hawaii.o $(CO2_AR)  
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/mpi_phase_space_integration_n2_ar.exe: examples/mpi_phase_space_integration_n2_ar.cpp $(MPI_OBJ) build/hep_hawaii.o $(N2_AR)  
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/mpi_phase_space_integration_ch4_co2.exe: examples/mpi_phase_space_integration_ch4_co2.cpp $(MPI_OBJ) build/hep_hawaii.o $(CH4_CO2)  
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/trajectory_co2_ar.exe: examples/trajectory_co2_ar.cpp build/trajectory.o $(OBJ) $(CO2_AR) 
	$(CXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/trajectory_h2_ar_requantized.exe: examples/trajectory_h2_ar_requantized.cpp build/trajectory.o $(OBJ) $(H2_AR)
	$(CXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/trajectory_ch4_co2.exe: examples/trajectory_ch4_co2.cpp build/trajectory.o $(OBJ) $(CH4_CO2) 
	$(CXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/correlation_co2_ar.exe: examples/correlation_co2_ar.cpp build/trajectory.o $(MPI_OBJ) $(CO2_AR) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) 

examples/correlation_n2_ar.exe: examples/correlation_n2_ar.cpp build/trajectory.o $(MPI_OBJ) $(N2_AR) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) 

examples/correlation_array_n2_ar.exe: examples/correlation_array_n2_ar.cpp build/trajectory.o $(MPI_OBJ) $(N2_AR)
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) 

examples/correlation_array_co2_ar.exe: examples/correlation_array_co2_ar.cpp build/trajectory.o $(MPI_OBJ) $(CO2_AR)
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) 

examples/prmu_calculation_co2_ar.exe: examples/prmu_calculation_co2_ar.cpp build/trajectory.o $(MPI_OBJ) $(CO2_AR) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) 

examples/fftrump.exe: examples/fftrump.cpp build/loess.o $(OBJ) 
	$(CXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm -lstdc++ $(LIB_SUNDIALS) $(LIB_GSL) -fopenmp

examples/correlation_co_ar.exe: examples/correlation_co_ar.cpp build/trajectory.o $(MPI_OBJ) $(CO_AR) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++ -lgfortran 
			
examples/prmu_calculation_co_ar_requantized.exe: examples/prmu_calculation_co_ar_requantized.cpp build/trajectory.o $(MPI_OBJ) $(CO_AR) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++ -lgfortran 

examples/prmu_calculation_h2_ar_requantized.exe: examples/prmu_calculation_h2_ar_requantized.cpp build/trajectory.o build/hep_hawaii.o $(MPI_OBJ) $(H2_AR) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++ -lgfortran 

examples/correlation_ch4_co2.exe: examples/correlation_ch4_co2.cpp build/trajectory.o $(MPI_OBJ) $(CH4_CO2) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) 

examples/correlation_array_ch4_co2.exe: examples/correlation_array_ch4_co2.cpp build/trajectory.o $(MPI_OBJ) $(CH4_CO2)
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) 

examples/correlation_he_ar.exe: examples/correlation_he_ar.cpp build/trajectory.o $(MPI_OBJ) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++ -lgfortran 

examples/test_sb.exe: examples/test_sb.c build/hawaii.o build/mtwist.o build/array.o build/trajectory.o
	$(CC) $(FLAGS) $(INC) -I./ $^ -o $@ -lm $(LIB_GSL) $(LIB_SUNDIALS)

examples/test_loess.exe: examples/test_loess.cpp build/hawaii.o build/mtwist.o build/array.o build/trajectory.o build/loess.o
	$(CC) $(FLAGS) $(INC) -fopenmp -I./ $^ -o $@ -lm $(LIB_GSL) $(LIB_SUNDIALS) -lstdc++

examples/test_fft.exe: examples/test_fft.c build/hawaii.o build/mtwist.o build/array.o build/trajectory.o build/loess.o
	$(CC) $(FLAGS) $(INC) -fopenmp -I./ $^ -o $@ -lm $(LIB_GSL) $(LIB_SUNDIALS) -lstdc++


build:
	mkdir -p $@

include $(wildcard build/*.d)
