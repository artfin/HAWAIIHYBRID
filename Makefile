.PHONY: all clean test docs

CC     ?= gcc
F      ?= gfortran
IFORT  ?= ifx
CXX    ?= g++
MPICC  ?= mpicc
MPICXX ?= mpic++

# set LD_PRELOAD before running program using sanitizer
# LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libasan.so.6 

# -Wswitch-enum: if default statement is present in the switch case, but not all the enum values are covered, the warning will still be emitted 
FLAGS_DEBUG   := -Wall -Wextra -Wswitch-enum -ggdb -O0 
FLAGS_RELEASE := -Wall -Wextra -Wswitch-enum -O2 -march=native -mtune=native # -pg -ggdb
#FLAGS_DEBUG   := -Wall -fsanitize=address -Wextra -Wswitch-enum -ggdb -O0 
#FLAGS_RELEASE := -Wall -fsanitize=address -Wextra -Wswitch-enum -O2 -march=native -mtune=native # -pg -ggdb
FLAGS_EIGEN   := -Wall -Wextra -Wswitch-enum -ggdb -O2 # -pg 
FLAGS := $(FLAGS_RELEASE)

LIB_GSL ?= -lgsl -lgslcblas 

-include Makefile.config

INC := $(INC_SUNDIALS) $(INC_EIGEN) $(INC_HEP) $(INC_GSL) -I./thirdparty/ -I./src/
LIBS := $(LIB_SUNDIALS) $(LIB_GSL) -lm -lstdc++ 

EXAMPLES := examples/phase_space_integration_co2_ar.exe      \
            examples/mpi_phase_space_integration_co2_ar.exe  \
            examples/mpi_phase_space_integration_ch4_co2.exe \
            examples/phase_space_integration_he_ar.exe       \
            examples/trajectory_co2_ar.exe                   \
            examples/trajectory_h2_ar_requantized.exe        \
            examples/trajectory_ch4_co2.exe                  \
            examples/trajectory_h2o_ar.exe                   \
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
            examples/prmu_calculation_line_test.exe          \
            examples/prmu_calculation_h2_ar_requantized.exe  \
            examples/prmu_calculation_d2_ar_requantized.exe  \
            examples/fftrump.exe                             \
            examples/test_sb.exe                             \
            examples/test_loess.exe                          \
            examples/test_fft.exe                            \
            driver.exe

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

build/hawaii.o: src/hawaii.c | build
	$(CC) $(FLAGS) $(INC) -c -MD $< -o $@

build/mpi_hawaii.o: src/hawaii.c | build
	$(MPICC) $(FLAGS) $(INC) -DUSE_MPI -c -MD $< -o $@

build/hep_hawaii.o: src/hep_hawaii.cpp | build
	$(MPICXX) $(FLAGS) $(INC) -c -MD $< -o $@

build/array.o: src/array.c | build
	$(CC) -std=c99 $(FLAGS) $(INC) -c -MD $< -o $@ 

build/trajectory.o: src/trajectory.c | build
	$(CC) -std=c99 $(FLAGS) $(INC) -c -MD $< -o $@ 

build/mtwist.o: thirdparty/mtwist.c | build
	$(CC) $(FLAGS) -c -MD $< -o $@ 

build/angles_handler.o: src/angles_handler.cpp | build
	$(CXX) $(FLAGS_EIGEN) $(INC) -c -MD -fPIC $< -o $@

build/loess.o: src/loess.cpp | build
	$(CXX) $(FLAGS_EIGEN) $(INC) -c -MD $< -o $@ 


###########################################################
####################### He-Ar #############################
###########################################################
build/HeAr.o: ./PES-IDS/HeAr.h | build
	$(CC) -DHEAR_IMPLEMENTATION -o $@ -x c -c $^ -lm -lstdc++

build/ai_pes_ids_he_ar.so: ./build/HeAr.o | build
	$(CC) -shared -o $@ $^ -lm -lstdc++

###########################################################
###################### N2-Ar-ISO #############################
###########################################################
build/n2_ar_pot_iso.o: ./PES-IDS/n2_ar_pot_iso.cpp | build
	$(CXX) $(FLAGS) $(INC) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 

build/n2_ar_pot_iso_der.o: ./PES-IDS/n2_ar_pot_iso_der.cpp | build
	$(CXX) $(FLAGS) $(INC) -c -MD -I./ $< -o $@ $(LINK_GSL) -lm 

###########################################################
###################### CO2-Ar #############################
###########################################################
build/ai_pes_co2_ar.o: ./PES-IDS/ai_pes_co2ar.c | build
	$(CC) $(FLAGS) $(INC) -c -MD -fPIC -I ./ $< -o $@ $(LINK_GSL) -lm 

build/ai_ids_co2_ar.o: ./PES-IDS/ai_ids_co2ar.cpp | build
	$(CXX) $(FLAGS) $(INC) -c -MD -fPIC -I ./ $< -o $@ $(LINK_GSL) -lm

build/ai_pes_co2ar_lib.o: ./PES-IDS/ai_pes_co2ar_lib.cpp
	$(CXX) $(INC) -c -MD -fPIC -I./ $< -o $@ -lm

build/ai_pes_co2ar.so: ./build/ai_pes_co2ar_lib.o build/ai_pes_co2_ar.o build/angles_handler.o | build
	$(CC) -shared -o $@ $^ -lm -lstdc++

build/ai_ids_co2ar_lib.o: ./PES-IDS/ai_ids_co2ar_lib.cpp
	$(CXX) $(INC) -c -MD -fPIC -I./ $< -o $@ -lm

build/ai_ids_co2ar.so: ./build/ai_ids_co2ar_lib.o build/ai_ids_co2_ar.o build/angles_handler.o | build
	$(CC) -shared -o $@ $^ -lm -lstdc++

###########################################################
###################### H2-Ar ##############################
###########################################################
build/c_basis_2_2_1_3_intermolecular.o: ./PES-IDS/c_basis_2_2_1_3_intermolecular.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm

build/c_basis_2_1_1_1_3_intermolecular.o: ./PES-IDS/c_basis_2_1_1_1_3_intermolecular.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm

build/c_basis_1_1_2_1_3_intermolecular.o: ./PES-IDS/c_basis_1_1_2_1_3_intermolecular.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm

build/ai_pes_h2ar_leroy.o: ./PES-IDS/ai_pes_h2ar_leroy.c | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ $(LINK_GSL) -lm

build/ai_ids_h2_ar_pip_nn.o: ./PES-IDS/ai_ids_h2_ar_pip_nn.cpp | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm

build/ai_pes_h2ar_leroy_lib.o: ./PES-IDS/ai_pes_h2ar_leroy_lib.cpp | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm 

build/ai_ids_h2_ar_pip_nn_lib.o: ./PES-IDS/ai_ids_h2_ar_pip_nn_lib.cpp | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm 

build/ai_ids_h2_ar_pip_nn.so: build/ai_ids_h2_ar_pip_nn_lib.o build/ai_ids_h2_ar_pip_nn.o \
							  build/c_basis_2_2_1_3_intermolecular.o build/c_basis_2_1_1_1_3_intermolecular.o build/c_basis_1_1_2_1_3_intermolecular.o \
							  build/angles_handler.o build/cnpy.o
	$(CC) -shared -o $@ $^ -lm -lstdc++

build/ai_pes_h2_ar_leroy.so: build/ai_pes_h2ar_leroy_lib.o build/ai_pes_h2ar_leroy.o build/angles_handler.o
	$(CC) -shared -o $@ $^ -lm

###########################################################
###################### CO-Ar ##############################
###########################################################
build/potv.o: ./PES-IDS/potv.f | build
	$(F) -c -fPIC $< -o $@ 

build/potv_d.o: ./PES-IDS/potv_d.f03 | build
	$(F) -c -fPIC $< -o $@

build/potv.so: ./PES-IDS/potv.cpp build/potv.o build/potv_d.o build/angles_handler.o
	$(CXX) $(INC_EIGEN) -shared -fPIC -I ./ -I./src/ -o $@ $^

build/ind_dipole_coar.so: ./PES-IDS/dipole_coar.cpp build/angles_handler.o 
	$(CXX) $(FLAGS) -shared -DDIPOLE_COAR_IMPLEMENTATION -I ./ -I ./src/ -fPIC $(INC_EIGEN) -o $@ $^ -lm

build/perm_dipole_coar.so: ./PES-IDS/perm_dipole_coar.c | build
	$(CC) $(FLAGS) -shared -I./ -I./src/ -o $@ $^ -lm

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
	$(CC) $(FLAGS) -c -MD -fPIC $< -o $@ -lm

build/ai_pes_n2_ar_pip_nn.o: ./PES-IDS/ai_pes_n2_ar_pip_nn.cpp | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c -MD -I./ $< -o $@ -lm

build/ai_ids_n2_ar_pip_nn.o: ./PES-IDS/ai_ids_n2_ar_pip_nn.cpp | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c -MD -I./ $< -o $@ -lm

###########################################################
##################### H2O-H2O #############################
###########################################################
build/c_basis_1_2_1_2_5_intermolecular.o: ./PES-IDS/c_basis_1_2_1_2_5_intermolecular.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm

build/c_jac_1_2_1_2_5_intermolecular.o: ./PES-IDS/c_jac_1_2_1_2_5_intermolecular.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm

build/ai_pes_h2o_h2o_nn_lib.o: ./PES-IDS/ai_pes_h2o_h2o_nn_lib.cpp | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm

build/ai_pes_h2o_h2o_nn.so: build/ai_pes_h2o_h2o_nn_lib.o \
							build/c_basis_1_2_1_2_5_intermolecular.o build/c_jac_1_2_1_2_5_intermolecular.o \
							build/angles_handler.o build/cnpy.o
	$(CC) -shared -o $@ $^ -lm -lstdc++ -lz

###########################################################
###################### H2O-Ar #############################
###########################################################
build/c_basis_1_2_1_4_purify.o: ./PES-IDS/c_basis_1_2_1_4_purify.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm

build/c_jac_1_2_1_4_purify.o: ./PES-IDS/c_jac_1_2_1_4_purify.cc | build
	$(CC) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm

build/ai_pes_h2o_ar_nn_lib.o: ./PES-IDS/ai_pes_h2o_ar_nn_lib.cpp | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm

build/ai_pes_h2o_ar_nn.so: build/ai_pes_h2o_ar_nn_lib.o \
							build/c_basis_1_2_1_4_purify.o build/c_jac_1_2_1_4_purify.o \
							build/angles_handler.o build/cnpy.o
	$(CC) -shared -o $@ $^ -lm -lstdc++ -lz

build/ai_ids_h2o_h2o_lib.o: ./PES-IDS/ai_ids_h2o_h2o_lib.cpp | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm

# H2O-H2O DMS Fortran sources (h4o2.dms4)
# Module dependency order: inv_share -> inv_mg321, inv_mg411 -> getdvec
# -r8 promotes default real to 8 bytes, matching the double precision
# arrays passed from calcdip (h4o2.dms4.f90) into getdvec and friends.
# -module keeps .mod files in build/ and searches there for them; without it
# they are written to and picked up from the working directory, which lets a
# stale .mod silently supply wrong mg321_ivs/mg321_ivb parameter values to
# getdvec (they are inlined at compile time).
H4O2_FFLAGS := -c -fPIC -r8 -module build -I./PES-IDS/

build/inv_share.o: ./PES-IDS/h2o-h2o/inv_share.f90 | build
	$(IFORT) $(H4O2_FFLAGS) $< -o $@

build/inv_mg321.o: ./PES-IDS/h2o-h2o/inv_mg321.f90 build/inv_share.o | build
	$(IFORT) $(H4O2_FFLAGS) $< -o $@

build/inv_mg411.o: ./PES-IDS/h2o-h2o/inv_mg411.f90 build/inv_share.o | build
	$(IFORT) $(H4O2_FFLAGS) $< -o $@

build/getdvec.o: ./PES-IDS/h2o-h2o/getdvec.f90 build/inv_share.o build/inv_mg321.o build/inv_mg411.o | build
	$(IFORT) $(H4O2_FFLAGS) $< -o $@

build/getd0.o: ./PES-IDS/h2o-h2o/getd0.f90 | build
	$(IFORT) $(H4O2_FFLAGS) $< -o $@

build/getr0.o: ./PES-IDS/h2o-h2o/getr0.f90 | build
	$(IFORT) $(H4O2_FFLAGS) $< -o $@

build/h4o2.dms4.o: ./PES-IDS/h2o-h2o/h4o2.dms4.f90 | build
	$(IFORT) $(H4O2_FFLAGS) $< -o $@

build/mgx_mk1d.o: ./PES-IDS/h2o-h2o/mgx_mk1d.f90 | build
	$(IFORT) $(H4O2_FFLAGS) $< -o $@

build/mgx_mk2d.o: ./PES-IDS/h2o-h2o/mgx_mk2d.f90 | build
	$(IFORT) $(H4O2_FFLAGS) $< -o $@

H4O2_DMS := build/inv_share.o build/inv_mg321.o build/inv_mg411.o \
            build/getdvec.o build/getd0.o build/getr0.o \
            build/h4o2.dms4.o build/mgx_mk1d.o build/mgx_mk2d.o

# -Bsymbolic-functions is load-bearing, not an optimization. ifx passes a
# contained procedure its host frame in %r10 (the SysV static-chain register),
# but exports the contained procedure globally, so the call to it goes through
# the PLT. glibc's lazy PLT resolver preserves the argument registers and
# clobbers %r10, so the first call to e.g. mg321_secs -> mg321_setd lands with a
# garbage host frame and segfaults. Binding intra-library calls directly removes
# the PLT indirection; -z now additionally resolves everything at load time.
build/ai_ids_h2o_h2o.so: build/ai_ids_h2o_h2o_lib.o build/angles_handler.o $(H4O2_DMS)
	$(IFORT) -shared -Wl,-Bsymbolic-functions -Wl,-z,now -o $@ $^ -lm -lstdc++

build/ai_ids_h2o_ar_nn_lib.o: ./PES-IDS/ai_ids_h2o_ar_nn_lib.cpp | build
	$(CXX) $(FLAGS) $(INC_EIGEN) -c -MD -fPIC -I./ $< -o $@ -lm

# H2O-Ar NN DMS (dipx/dipy/dipz). Fixed-form F77, hence -std=legacy.
# Only dipx.f defines tranfun and includes dms_interface.f; dipy/dipz call tranfun
# as an external, so dipx.o must always be linked in alongside them.
# The nets read their weights from PES-IDS/h2o-ar-dms/ using a path relative to the
# working directory, so the driver has to be run from the repository root.
H2OAR_DMS_FFLAGS := -c -fPIC -std=legacy -O2 -I./PES-IDS/h2o-ar-dms/ -J build

build/dipx.o: ./PES-IDS/h2o-ar-dms/dipx.f ./PES-IDS/h2o-ar-dms/dms_interface.f | build
	$(F) $(H2OAR_DMS_FFLAGS) $< -o $@

build/dipy.o: ./PES-IDS/h2o-ar-dms/dipy.f | build
	$(F) $(H2OAR_DMS_FFLAGS) $< -o $@

build/dipz.o: ./PES-IDS/h2o-ar-dms/dipz.f | build
	$(F) $(H2OAR_DMS_FFLAGS) $< -o $@

H2OAR_DMS := build/dipx.o build/dipy.o build/dipz.o

build/ai_ids_h2o_ar_nn.so: build/ai_ids_h2o_ar_nn_lib.o build/angles_handler.o $(H2OAR_DMS)
	$(CC) -shared -o $@ $^ -lm -lstdc++ -lgfortran

###########################################################
##################### CH4-CO2 #############################
###########################################################
build/ai_pes_ch4_co2.o: ./PES-IDS/ai_pes_ch4_co2.c | build
	$(CC) $(FLAGS) $(INC_GSL) -c -MD -fPIC -I./ $< -o $@ $(LINK_GSL) -lm 

build/ai_pes_ch4_co2_dEdR.o: ./PES-IDS/ai_pes_ch4_co2_dEdR.c | build
	$(CC) $(FLAGS) $(INC_GSL) -c -MD -fPIC -I./ $< -o $@ $(LINK_GSL) -lm 

build/ai_pes_ch4_co2_dEdphi1.o: ./PES-IDS/ai_pes_ch4_co2_dEdphi1.c | build
	$(CC) $(FLAGS) $(INC_GSL) -c -MD -fPIC -I./ $< -o $@ $(LINK_GSL) -lm 

build/ai_pes_ch4_co2_dEdtheta1.o: ./PES-IDS/ai_pes_ch4_co2_dEdtheta1.c | build
	$(CC) $(FLAGS) $(INC_GSL) -c -MD -fPIC -I./ $< -o $@ $(LINK_GSL) -lm 

build/ai_pes_ch4_co2_dEdphi2.o: ./PES-IDS/ai_pes_ch4_co2_dEdphi2.c | build
	$(CC) $(FLAGS) $(INC_GSL) -c -MD -fPIC -I./ $< -o $@ $(LINK_GSL) -lm 

build/ai_pes_ch4_co2_dEdtheta2.o: ./PES-IDS/ai_pes_ch4_co2_dEdtheta2.c | build
	$(CC) $(FLAGS) $(INC_GSL) -c -MD -fPIC -I./ $< -o $@ $(LINK_GSL) -lm 

build/ai_pes_ch4_co2_lib.o: ./PES-IDS/ai_pes_ch4_co2_lib.cpp | build
	$(CXX) $(FLAGS) $(INC) -c -MD -fPIC -I./ -I./PES-IDS/ $< -o $@ -lm

build/ai_pes_ch4_co2.so: ./build/ai_pes_ch4_co2_lib.o build/angles_handler.o \
	  					 ./build/ai_pes_ch4_co2.o ./build/ai_pes_ch4_co2_dEdR.o \
						 ./build/ai_pes_ch4_co2_dEdphi1.o build/ai_pes_ch4_co2_dEdtheta1.o \
						 ./build/ai_pes_ch4_co2_dEdphi2.o build/ai_pes_ch4_co2_dEdtheta2.o 
	$(CC) -shared -o $@ $^ -lm

build/ai_ids_ch4_co2.o: ./PES-IDS/ai_ids_ch4_co2.cpp | build
	$(CC) $(FLAGS) $(INC_GSL) $(INC_EIGEN) -c -fPIC -MD -I./ $< -o $@ $(LINK_GSL) -lm 

build/ai_ids_ch4_co2_lib.o: ./PES-IDS/ai_ids_ch4_co2_lib.cpp | build
	$(CC) $(FLAGS) $(INC_GSL) $(INC_EIGEN) -c -fPIC -MD -I./ $< -o $@ $(LINK_GSL) -lm 

build/ai_ids_ch4_co2.so: ./build/ai_ids_ch4_co2.o ./build/ai_ids_ch4_co2_lib.o ./build/angles_handler.o
	$(CC) -shared -o $@ $^ -lm
###########################################################

OBJ     := build/hawaii.o build/mtwist.o build/angles_handler.o build/array.o build/trajectory.o
MPI_OBJ := build/mpi_hawaii.o build/mtwist.o build/angles_handler.o build/array.o build/trajectory.o build/hep_hawaii.o
CO2_AR  := build/ai_pes_co2_ar.o build/ai_ids_co2_ar.o build/ai_pes_co2ar_lib.o build/ai_ids_co2ar_lib.o
N2_AR_ISO  := build/n2_ar_pot_iso.o build/n2_ar_pot_iso_der.o build/cnpy.o -lz  build/ai_ids_n2_ar_pip_nn.o \
		   build/c_basis_2_2_1_3_purify.o build/c_basis_2_1_1_1_3_purify.o build/c_basis_1_1_2_1_3_purify.o
N2_AR   := build/cnpy.o -lz build/ai_pes_n2_ar_pip_nn.o build/ai_ids_n2_ar_pip_nn.o \
		   build/c_basis_2_1_4_purify.o build/c_jac_2_1_4_purify.o \
		   build/c_basis_2_2_1_3_purify.o build/c_basis_2_1_1_1_3_purify.o build/c_basis_1_1_2_1_3_purify.o
H2_AR   := build/cnpy.o -lz build/ai_pes_h2ar_leroy_lib.o build/ai_pes_h2ar_leroy.o build/ai_ids_h2_ar_pip_nn.o \
		   build/c_basis_2_2_1_3_intermolecular.o build/c_basis_2_1_1_1_3_intermolecular.o build/c_basis_1_1_2_1_3_intermolecular.o
CO_AR   := build/potv.o build/potv_d.o 
H2O_AR  := build/cnpy.o -lz build/ai_pes_h2o_ar_nn_lib.o \
		   build/c_basis_1_2_1_4_purify.o build/c_jac_1_2_1_4_purify.o
CH4_CO2 := build/ai_pes_ch4_co2_lib.o build/ai_ids_ch4_co2_lib.o \
		   build/ai_pes_ch4_co2.o build/ai_pes_ch4_co2_dEdR.o \
		   build/ai_pes_ch4_co2_dEdphi1.o build/ai_pes_ch4_co2_dEdtheta1.o \
		   build/ai_pes_ch4_co2_dEdphi2.o build/ai_pes_ch4_co2_dEdtheta2.o \
		   build/ai_ids_ch4_co2.o
 
 
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

examples/trajectory_h2o_ar.exe: examples/trajectory_h2o_ar.cpp build/trajectory.o $(OBJ) $(H2O_AR)
	$(CXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/trajectory_ch4_co2.exe: examples/trajectory_ch4_co2.cpp build/trajectory.o $(OBJ) $(CH4_CO2) 
	$(CXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++  

examples/correlation_co2_ar.exe: examples/correlation_co2_ar.cpp build/trajectory.o $(MPI_OBJ) $(CO2_AR) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) 

examples/correlation_n2_ar_iso.exe: examples/correlation_n2_ar_iso.cpp build/trajectory.o $(MPI_OBJ) $(N2_AR_ISO) 
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
	$(CXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm -lstdc++ $(LIB_SUNDIALS) $(LIB_GSL) 

examples/correlation_co_ar.exe: examples/correlation_co_ar.cpp build/trajectory.o $(MPI_OBJ) $(CO_AR) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++ -lgfortran 
			
examples/prmu_calculation_co_ar_requantized.exe: examples/prmu_calculation_co_ar_requantized.cpp build/trajectory.o $(MPI_OBJ) $(CO_AR) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++ -lgfortran 

examples/prmu_calculation_h2_ar_requantized.exe: examples/prmu_calculation_h2_ar_requantized.cpp build/trajectory.o $(MPI_OBJ) $(H2_AR) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++ -lgfortran 

examples/prmu_calculation_d2_ar_requantized.exe: examples/prmu_calculation_d2_ar_requantized.cpp build/trajectory.o $(MPI_OBJ) $(H2_AR) 
	$(MPICXX) $(FLAGS) $(INC) -I./ -I./PES-IDS/ $^ -o $@ -lm $(LIB_SUNDIALS) $(LIB_GSL) -lstdc++ -lgfortran 

examples/prmu_calculation_line_test.exe: examples/prmu_calculation_line_test.cpp build/trajectory.o $(MPI_OBJ)
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
	$(CC) $(FLAGS) $(INC) -I./ $^ -o $@ -lm $(LIB_GSL) $(LIB_SUNDIALS) -lstdc++

examples/test_fft.exe: examples/test_fft.c build/hawaii.o build/mtwist.o build/array.o build/trajectory.o build/loess.o
	$(CC) $(FLAGS) $(INC) -I./ $^ -o $@ -lm $(LIB_GSL) $(LIB_SUNDIALS) -lstdc++

GIT_COMMIT := $(shell git rev-parse --short HEAD 2>/dev/null || echo "unknown")
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")

driver.exe: src/driver.c build/mpi_hawaii.o build/mtwist.o build/trajectory.o build/array.o build/angles_handler.o build/hep_hawaii.o build/loess.o
	$(MPICC) -Wall -Wextra -ggdb $(INC) -DGIT_COMMIT='"$(GIT_COMMIT)"' -DGIT_BRANCH='"$(GIT_BRANCH)"' $^ -o $@ -lm $(LIB_GSL) $(LIB_SUNDIALS) -lstdc++ -ldl -lpthread # -lasan 

hawaii_test.exe: src/hawaii_test.c driver.exe
	$(CC) $(FLAGS) -isystem ./thirdparty/ $< -o $@

d4b.exe: d4b.c $(OBJ)
	$(CC) $(FLAGS) $(INC) $^ -o $@ $(LIBS) #-lasan 


build:
	mkdir -p $@

include $(wildcard build/*.d)
