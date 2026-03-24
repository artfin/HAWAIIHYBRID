// Shared library wrapper for H2O-H2O induced dipole surface
// Exposes dipole_lab, dipole_init for HAWAIIHYBRID
//
// Coordinate mapping:
//   HAWAIIHYBRID intermediate_q = [Phi, Theta, R, phi1T, theta1T, psi1T, phi2T, theta2T, psi2T]
//   PES code r_ang              = [R,   Phi, Theta, phi1T, theta1T, psi1T, phi2T, theta2T, psi2T]

#include <cstdio>

#include "src/angles_handler.hpp"

// Fortran dipole surface (libh4o2dms4ifc.a)
extern "C" {
    void predip_();
    void calcdip_(double *dip, double *cart_in);
}

extern "C" {

void dipole_init(bool log)
{
    predip_();
    if (log) {
        printf("H2O-H2O dipole surface initialized\n");
    }
}

void dipole_lab(double *q, double diplab[3])
{
    double r_ang[9] = {q[2], q[0], q[1], q[3], q[4], q[5], q[6], q[7], q[8]};

    double cart[3][6];
    h2o_h2o_lab_to_cart(r_ang, cart);

    calcdip_(diplab, &cart[0][0]);
}

} // extern "C"
