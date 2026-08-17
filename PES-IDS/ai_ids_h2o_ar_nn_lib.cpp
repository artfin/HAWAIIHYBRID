// Shared-library wrapper for the FULL H2O-Ar NN dipole moment surface.
// Exposes dipole_lab and dipole_init for HAWAIIHYBRID.
//
// Model: the three-component NN in PES-IDS/h2o-ar-dms (dipx.f, dipy.f, dipz.f),
// from Liu, Wang, Zhou, and Xie, "A Full-Dimensional ab initio Intermolecular
// Potential Energy Surface and Dipole Moment Surfaces for H2O-Ar", Current
// Chinese Science 2, 325 (2022), doi:10.2174/2210298102666220404103308.
//
// IMPORTANT SEMANTICS: the fitted surface is a FULL/TOTAL dipole surface.  It
// includes the permanent dipole of the H2O monomer and approaches that nonzero
// monomer contribution as Ar is separated.  Accordingly, this wrapper returns
//
//     mu_full = mu_Liu,H2O-Ar
//
// without subtracting any monomer reference.  It must not be called an induced
// dipole surface.  For the rigid-monomer induced quantity use the separate target
// build/ind_dipole_h2o_ar.so, which returns
//
//     mu_ind(R,Omega) = mu_full(R,Omega) - mu_full(20 bohr,Omega).
//
// The second term is the monomer contribution embedded in the Liu fit, evaluated
// at the upper edge of its R training range with the same angular coordinates.
// The standalone build/perm_dipole_h2o_ar.so returns this identical reference,
// allowing total, induced, permanent, and mixed contributions to be selected in
// the same manner as for CO-Ar.
//
// Electronic-structure/fitting details: the Liu DMS is MP2/AVTZ and was fitted
// over R=4-20 bohr with flexible Radau H2O coordinates.  Each Cartesian component
// is an average of three 9->50->50->1 neural networks.  Output is in atomic units.
// The raw networks use frame-dependent Cartesian inputs; their fitted-frame
// reconstruction and hydrogen-labelling logic live in h2o_ar_nn_dms_common.cpp.
// No damping or long-range correction is applied outside the fitted R interval.
//
// Coordinate mapping:
//   HAWAIIHYBRID q = [Phi, Theta, R, phi1T, theta1T, psi1T]
//   PES r_ang      = [R,   Phi,   Theta, phi1T, theta1T, psi1T]
// Monomer 1 is rigid H2O and monomer 2 is Ar.

#include <cstdio>

#include "PES-IDS/h2o_ar_nn_dms_common.hpp"

extern "C" {

void dipole_init(bool log)
{
    h2o_ar_nn_dms::initialize();
    if (log) {
        std::printf("H2O-Ar full Liu NN dipole surface initialized "
                    "(includes permanent H2O dipole)\n");
    }
}

void dipole_lab(double *q, double diplab[3])
{
    h2o_ar_nn_dms::total_dipole_lab(q, diplab);
}

} // extern "C"
