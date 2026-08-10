// Shared library wrapper for the H2O-H2O "fi-nn" neural network PES
// (vendored from water/2026/fi_nn_smoke_test/cpp_export_fast).
// Exposes pes_lab, dpes_lab, pes_init for HAWAIIHYBRID
//
// Coordinate mapping:
//   HAWAIIHYBRID intermediate_q = [Phi, Theta, R, phi1T, theta1T, psi1T, phi2T, theta2T, psi2T]
//   PES code r_ang              = [R,   Phi, Theta, phi1T, theta1T, psi1T, phi2T, theta2T, psi2T]
//
// Atom order mismatch: h2o_h2o_lab_to_cart() produces cart[3][6] in this
// repo's order (Oa, Ha1, Ha2, Ob, Hb1, Hb2); fi_nn_pes expects Cartesian
// input in "Dalian" order (Oa, Ob, Ha1, Ha2, Hb1, Hb2). Both use atom-major
// flattening (index = 3*atom + dim), so bridging the two is a fixed
// per-atom-block permutation -- FINN_ATOM_FROM_REPO[d] gives the repo atom
// index that fills Dalian slot d.

#include <Eigen/Dense>

#include "src/constants.h"
#include "src/angles_handler.hpp"

#include "fi_nn_pes.hpp"

#define MODEL_NPZ "PES-IDS/npz/water-finn.npz"

static const int FINN_ATOM_FROM_REPO[6] = {0, 3, 1, 2, 4, 5};

// Short-range repulsive-wall guard: DISABLED for now.
//
// This was copied verbatim from the pip-nn H2O-H2O PES (ai_pes_h2o_h2o_nn_lib.cpp),
// whose short-range guard was a palliative for a well-documented pip-nn failure
// mode (spurious deep wells, up to -56,000 cm-1, concentrated right around
// R~4.2 Bohr -- see water/2026/lab_notes.md Sec.19.10). fi-nn's own 60-random-
// orientation stress test (Sec.19.9) found nothing resembling that: worst-case
// spurious wobble 9.92 cm-1, sub-visual, vs. pip-nn's up to -56,000 cm-1.
// Independently, a 1D radial scan and a ~3900-trial multi-start search here
// (unguarded) found the raw fi-nn surface rises smoothly/monotonically all the
// way from the genuine global minimum (closest contact 3.75 Bohr, -1659 cm-1)
// through the physical repulsive wall out to atom-overlap distances, with zero
// spurious close-contact wells detected.
//
// That said, none of this is a rigorous training-data-coverage check (which
// would be the principled way to set a cutoff), so the guard is kept here,
// commented out rather than deleted, so it can be reinstated (with a properly
// justified threshold) if evidence of bad short-range behavior ever turns up.
//
// static const double RMIN_INTER_SQ = 4.0 * 4.0; // 4 Bohr ~ 2.1 A
//
// static bool inter_monomer_too_close(double cart[3][6])
// {
//     // Monomer A: atoms 0,1,2  Monomer B: atoms 3,4,5
//     for (int a = 0; a < 3; ++a) {
//         for (int b = 3; b < 6; ++b) {
//             double dx = cart[0][a] - cart[0][b];
//             double dy = cart[1][a] - cart[1][b];
//             double dz = cart[2][a] - cart[2][b];
//             if (dx*dx + dy*dy + dz*dz < RMIN_INTER_SQ)
//                 return true;
//         }
//     }
//     return false;
// }

static void cart_to_finn_cart18(double cart[3][6], double cart18[18])
{
    for (int d = 0; d < 6; ++d) {
        int r = FINN_ATOM_FROM_REPO[d];
        cart18[3*d + 0] = cart[0][r];
        cart18[3*d + 1] = cart[1][r];
        cart18[3*d + 2] = cart[2][r];
    }
}

extern "C" {

void pes_init(void)
{
    fi_nn_pes_init(MODEL_NPZ);
}

double pes_lab(double *q)
{
    double r_ang[9] = {q[2], q[0], q[1], q[3], q[4], q[5], q[6], q[7], q[8]};

    // Short-range guard disabled -- see comment above inter_monomer_too_close().
    // if (r_ang[0] <= 4.0)
    //     return 1e5;

    double cart[3][6];
    h2o_h2o_lab_to_cart(r_ang, cart);

    // if (inter_monomer_too_close(cart))
    //     return 1e5;

    double cart18[18];
    cart_to_finn_cart18(cart, cart18);

    // fi_nn_pes_energy returns energy in eV, convert to Hartree
    return fi_nn_pes_energy(cart18) / HTOEV;
}

void dpes_lab(double *q, double *dpesdq)
{
    double r_ang[9] = {q[2], q[0], q[1], q[3], q[4], q[5], q[6], q[7], q[8]};

    // Short-range guard disabled -- see comment above inter_monomer_too_close().
    // if (r_ang[0] <= 4.0) {
    //     for (int i = 0; i < 9; ++i)
    //         dpesdq[i] = 0.0;
    //     return;
    // }

    double cart[3][6];
    h2o_h2o_lab_to_cart(r_ang, cart);

    // if (inter_monomer_too_close(cart)) {
    //     for (int i = 0; i < 9; ++i)
    //         dpesdq[i] = 0.0;
    //     return;
    // }

    double cart18[18];
    cart_to_finn_cart18(cart, cart18);

    // Get Cartesian gradient from fi-nn (eV / Bohr, Dalian atom order)
    double energy_eV;
    double grad18[18];
    fi_nn_pes_energy_and_gradient(cart18, &energy_eV, grad18);

    // Un-permute back into repo atom order and convert to Hartree / Bohr
    Eigen::VectorXd dVdcart(18);
    for (int d = 0; d < 6; ++d) {
        int r = FINN_ATOM_FROM_REPO[d];
        for (int dim = 0; dim < 3; ++dim)
            dVdcart(3*r + dim) = grad18[3*d + dim] / HTOEV;
    }

    // Jacobian d(cart) / d(r_ang), shape 9x18
    Eigen::MatrixXd jacobi_mat(9, 18);
    jacobi_mat.setZero(9, 18);
    h2o_h2o_der_cart_by_rang(jacobi_mat, cart, r_ang);

    // dV/d(r_ang) = J * dV/d(cart)
    Eigen::Matrix<double, 9, 1> dVdrang = jacobi_mat * dVdcart;

    // Remap from r_ang order [dV/dR, dV/dPhi, dV/dTheta, ...]
    // to intermediate_q order [dV/dPhi, dV/dTheta, dV/dR, ...]
    dpesdq[0] = dVdrang(1); // dV/dPhi
    dpesdq[1] = dVdrang(2); // dV/dTheta
    dpesdq[2] = dVdrang(0); // dV/dR
    dpesdq[3] = dVdrang(3); // dV/dphi1T
    dpesdq[4] = dVdrang(4); // dV/dtheta1T
    dpesdq[5] = dVdrang(5); // dV/dpsi1T
    dpesdq[6] = dVdrang(6); // dV/dphi2T
    dpesdq[7] = dVdrang(7); // dV/dtheta2T
    dpesdq[8] = dVdrang(8); // dV/dpsi2T
}

} // extern "C"
