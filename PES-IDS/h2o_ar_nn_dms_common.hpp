#ifndef H2O_AR_NN_DMS_COMMON_HPP
#define H2O_AR_NN_DMS_COMMON_HPP

namespace h2o_ar_nn_dms {

// Upper edge of the Liu NN training interval.  At fixed rigid-H2O orientation,
// the full dipole evaluated here is used as the embedded monomer reference.
constexpr double monomer_reference_R_bohr = 20.0;

// Initialize the three Fortran neural-network components.  The caller must invoke
// this exactly once for each loaded shared-library instance before evaluation.
void initialize();

// Evaluate the full Liu et al. H2O-Ar dipole in the laboratory frame.  q follows
// the HAWAIIHYBRID ordering [Phi, Theta, R, phi1T, theta1T, psi1T].
void total_dipole_lab(double *q, double diplab[3]);

// Evaluate the monomer term embedded in the Liu surface by retaining all angular
// coordinates in q and replacing only R with monomer_reference_R_bohr.
void embedded_monomer_dipole_lab(double *q, double diplab[3]);

} // namespace h2o_ar_nn_dms

#endif // H2O_AR_NN_DMS_COMMON_HPP
