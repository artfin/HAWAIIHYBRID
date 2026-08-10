.PHONY: all clean

# Standalone build for the H2O-H2O dipole moment surface (Fortran DMS +
# thin C++ wrapper). This remains available in addition to the integrated root
# Makefile target. Run it directly (`make -C PES-IDS/h2o-h2o-dipole`), overriding IFORT/
# INC_EIGEN on the command line or via ../../Makefile.config as needed.
# Sources are referenced in place (../ai_ids_h2o_h2o_lib.cpp, ../h2o-h2o/
# *.f90) -- nothing here is a copy/fork of the production dipole code, and
# the output lands at the SAME conventional path
# (../../build/ai_ids_h2o_h2o.so) the root Makefile's own dipole target
# produces, so example configs (so_dipole = "build/ai_ids_h2o_h2o.so")
# don't need to know which Makefile built it.
#
# predip() loads PES-IDS/h2o-h2o/h4o2.dms4.coeff.dat at runtime, matching the
# repository's other data paths. Run consumers from the repository root.

CXX   ?= g++
IFORT ?= ifx

-include ../../Makefile.config

FLAGS := -Wall -Wextra -O2 -fPIC
INC   := $(INC_EIGEN) -I../..

TARGET := ../../build/ai_ids_h2o_h2o.so

# Module dependency order: inv_share -> inv_mg321, inv_mg411 -> getdvec
# (all objects land in this directory, so .mod files resolve via the
# default cwd search -- no -I needed for module lookup).
H4O2_FFLAGS := -c -fPIC -r8

DMS_OBJS := inv_share.o inv_mg321.o inv_mg411.o getdvec.o getd0.o getr0.o h4o2.dms4.o mgx_mk1d.o mgx_mk2d.o

all: $(TARGET)

../../build:
	mkdir -p ../../build

$(TARGET): ai_ids_h2o_h2o_lib.o angles_handler.o $(DMS_OBJS) | ../../build
	$(IFORT) -shared -o $@ $^ -lstdc++ -lm

ai_ids_h2o_h2o_lib.o: ../ai_ids_h2o_h2o_lib.cpp
	$(CXX) $(FLAGS) $(INC) -c $< -o $@

angles_handler.o: ../../src/angles_handler.cpp
	$(CXX) $(FLAGS) $(INC) -c $< -o $@

inv_share.o: ../h2o-h2o/inv_share.f90
	$(IFORT) $(H4O2_FFLAGS) $< -o $@
inv_mg321.o: ../h2o-h2o/inv_mg321.f90 inv_share.o
	$(IFORT) $(H4O2_FFLAGS) $< -o $@
inv_mg411.o: ../h2o-h2o/inv_mg411.f90 inv_share.o
	$(IFORT) $(H4O2_FFLAGS) $< -o $@
getdvec.o: ../h2o-h2o/getdvec.f90 inv_share.o inv_mg321.o inv_mg411.o
	$(IFORT) $(H4O2_FFLAGS) $< -o $@
getd0.o: ../h2o-h2o/getd0.f90
	$(IFORT) $(H4O2_FFLAGS) $< -o $@
getr0.o: ../h2o-h2o/getr0.f90
	$(IFORT) $(H4O2_FFLAGS) $< -o $@
h4o2.dms4.o: ../h2o-h2o/h4o2.dms4.f90 ../h2o-h2o/h4o2.dms4.coeff.dat
	$(IFORT) $(H4O2_FFLAGS) $< -o $@
mgx_mk1d.o: ../h2o-h2o/mgx_mk1d.f90
	$(IFORT) $(H4O2_FFLAGS) $< -o $@
mgx_mk2d.o: ../h2o-h2o/mgx_mk2d.f90
	$(IFORT) $(H4O2_FFLAGS) $< -o $@

clean:
	rm -f *.o *.mod $(TARGET)
