// Monomer dipole embedded in the full Liu H2O-Ar NN surface.
//
// This is the H2O-Ar analogue of perm_dipole_coar.c.  It exports the standard
// HAWAIIHYBRID dipole_lab ABI and evaluates the full Liu DMS at R_ref=20 bohr,
// the upper edge of its fitted interval.  The input H2O/angle coordinates are
// retained and the actual pair separation is ignored.
//
// Ideally this large-R result would be one rigid-H2O vector independent of the Ar
// direction.  In practice the Cartesian NN fit has a small orientation-dependent
// asymptotic error.  Returning the same fitted reference used by
// ind_dipole_h2o_ar.cpp guarantees full = induced + permanent component by component.

#include <cstdio>

#include "PES-IDS/h2o_ar_nn_dms_common.hpp"

extern "C" {

void dipole_init(bool log)
{
    h2o_ar_nn_dms::initialize();
    if (log) {
        std::printf("Embedded Liu H2O monomer dipole initialized "
                    "(full DMS at R=%.1f bohr)\n",
                    h2o_ar_nn_dms::monomer_reference_R_bohr);
    }
}

void dipole_lab(double *q, double diplab[3])
{
    h2o_ar_nn_dms::embedded_monomer_dipole_lab(q, diplab);
}

} // extern "C"
