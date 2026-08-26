// Permanent dipole of the isolated rigid H2O monomer.
//
// HAWAIIHYBRID coordinates are
//   q = [Phi, Theta, R, phi1T, theta1T, psi1T].
// Only the H2O Euler angles q[3..5] describe the monomer orientation.  The
// intermolecular coordinates q[0..2] are intentionally ignored, so neither the
// position nor the direction of Ar can affect this observable.
//
// At the rigid geometry used by the H2O-Ar trajectory code, direct evaluation
// of the Schwenke-Partridge 2000 DMS gives the body-fixed dipole
//
//   mu_body = (0, 0, -0.7359507568) a.u.
//
// The minus sign follows the convention in src/constants.h: H2O lies in the xz
// plane, its C2 axis is z, and the hydrogen side points toward body-fixed -z.

#include <cmath>
#include <cstdio>

namespace {

constexpr double mu_sp2000 = 0.7359507568;

} // namespace

extern "C" {

void dipole_init(bool log)
{
    if (log) {
        std::printf("Rigid H2O SP2000 permanent dipole initialized "
                    "(magnitude %.10f a.u.; Ar coordinates ignored)\n",
                    mu_sp2000);
    }
}

void dipole_lab(double *q, double diplab[3])
{
    // The H2O Cartesian geometry uses
    //   S1 = Sz(phi1T)^T Sx(theta1T)^T Sz(psi1T)^T.
    // Since mu_body lies on z, the final rotation about body-fixed z (psi1T)
    // leaves it unchanged.  Expanding S1 * mu_body gives the expressions below.
    const double phi1T = q[3];
    const double theta1T = q[4];
    const double sin_theta1T = std::sin(theta1T);

    diplab[0] = -mu_sp2000 * std::sin(phi1T) * sin_theta1T;
    diplab[1] =  mu_sp2000 * std::cos(phi1T) * sin_theta1T;
    diplab[2] = -mu_sp2000 * std::cos(theta1T);
}

} // extern "C"
