// H2O-Ar interaction-induced dipole wrapper for rigid H2O.
//
// Definition:
//   mu_ind(R,Omega) = mu_Liu,total(R,Omega) - mu_Liu,total(R_ref,Omega)
//
// Both terms come from the same three-component NN DMS of Liu et al. (2022).
// R_ref=20 bohr is the largest separation in the fitted interval and represents
// the monomer dipole embedded in that surface.  All H2O and angular coordinates
// are identical between the two evaluations; only R is replaced by R_ref.
//
// This wrapper is intentionally rigid-monomer-only: q contains no H2O vibrational
// coordinates, and both terms reconstruct H2O through h2o_ar_lab_to_cart.
// Subtracting the NN's own large-R value cancels its orientation-dependent error in
// the fitted monomer limit and makes mu_ind exactly zero at R=R_ref.  No evaluation
// beyond the trained R range and no artificial cutoff or damping are used.

#include <cstdio>

#include "PES-IDS/h2o_ar_nn_dms_common.hpp"

extern "C" {

void dipole_init(bool log)
{
    h2o_ar_nn_dms::initialize();
    if (log) {
        std::printf("H2O-Ar induced dipole initialized "
                    "(full Liu DMS minus its embedded monomer at R=%.1f bohr)\n",
                    h2o_ar_nn_dms::monomer_reference_R_bohr);
    }
}

void dipole_lab(double *q, double diplab[3])
{
    double total[3], monomer[3];
    h2o_ar_nn_dms::total_dipole_lab(q, total);
    h2o_ar_nn_dms::embedded_monomer_dipole_lab(q, monomer);

    for (int axis = 0; axis < 3; ++axis) {
        diplab[axis] = total[axis] - monomer[axis];
    }
}

} // extern "C"
