.PHONY: all clean

# Standalone build for the vendored "fi-nn" H2O-H2O PES.
# This remains available in addition to the integrated root Makefile target.
# Run it directly (`make -C PES-IDS/h2o-h2o-finn`) or from the repo root
# (`make -C PES-IDS/h2o-h2o-finn all`). Reuses the repo's own Eigen include
# path from ../../Makefile.config (see ../../Makefile.config.example) and its
# already-vendored cnpy.{h,cpp} / angles_handler.{hpp,cpp} -- nothing else
# from the root build is required to exist first.

CXX ?= g++

-include ../../Makefile.config

FLAGS := -Wall -Wextra -O2 -fPIC
INC   := $(INC_EIGEN) -I../..

TARGET := ../../build/ai_pes_h2o_h2o_finn.so
OBJS   := weights.o fi_nn_pes.o cnpy.o angles_handler.o ai_pes_h2o_h2o_finn_lib.o

all: $(TARGET)

../../build:
	mkdir -p ../../build

$(TARGET): $(OBJS) | ../../build
	$(CXX) -shared -fPIC -o $@ $(OBJS) -lm -lz

weights.o: weights.cpp weights.hpp ../cnpy.h
	$(CXX) $(FLAGS) $(INC) -c weights.cpp -o weights.o

fi_nn_pes.o: fi_nn_pes.cpp fi_nn_pes.hpp fast_pes_pipeline.hpp orbit_basis_generated.hpp features_jacobian_generated.hpp weights.hpp
	$(CXX) $(FLAGS) $(INC) -c fi_nn_pes.cpp -o fi_nn_pes.o

cnpy.o: ../cnpy.cpp ../cnpy.h
	$(CXX) $(FLAGS) $(INC) -c ../cnpy.cpp -o cnpy.o

angles_handler.o: ../../src/angles_handler.cpp ../../src/angles_handler.hpp
	$(CXX) $(FLAGS) $(INC) -c ../../src/angles_handler.cpp -o angles_handler.o

ai_pes_h2o_h2o_finn_lib.o: ai_pes_h2o_h2o_finn_lib.cpp fi_nn_pes.hpp ../../src/constants.h ../../src/angles_handler.hpp
	$(CXX) $(FLAGS) $(INC) -c ai_pes_h2o_h2o_finn_lib.cpp -o ai_pes_h2o_h2o_finn_lib.o

clean:
	rm -f *.o $(TARGET)
