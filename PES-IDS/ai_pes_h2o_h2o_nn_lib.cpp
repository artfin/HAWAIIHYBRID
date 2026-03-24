// Shared library wrapper for H2O-H2O neural network PES
// Exposes pes_lab, dpes_lab, pes_init for HAWAIIHYBRID
//
// Coordinate mapping:
//   HAWAIIHYBRID intermediate_q = [Phi, Theta, R, phi1T, theta1T, psi1T, phi2T, theta2T, psi2T]
//   PES code r_ang              = [R,   Phi, Theta, phi1T, theta1T, psi1T, phi2T, theta2T, psi2T]

#include <cstring>

#include <Eigen/Dense>

#include "src/constants.h"
#include "src/angles_handler.hpp"

#include "PES-IDS/c_basis_1_2_1_2_5_intermolecular.h"
#include "PES-IDS/c_jac_1_2_1_2_5_intermolecular.h"

static const int natoms = 6;
static const int ndist = natoms * (natoms - 1) / 2;
static const double a0 = 2.0;

#define MODEL_NPZ "PES-IDS/npz/water-extended.npz"

// ---- PIP coordinate transforms ----

static void make_yij_1_2_1_2_5_intermolecular(const double *x, double *yij, int natoms)
{
    double drx, dry, drz;

    size_t k = 0;
    for (size_t i = 0; i < (size_t)natoms; ++i) {
        for (size_t j = i + 1; j < (size_t)natoms; ++j) {
            if (i == 0 && j == 1) { yij[k] = 0.0; k = k + 1; continue; } // O H1
            if (i == 0 && j == 2) { yij[k] = 0.0; k = k + 1; continue; } // O H2
            if (i == 1 && j == 2) { yij[k] = 0.0; k = k + 1; continue; } // H1 H2
            if (i == 3 && j == 4) { yij[k] = 0.0; k = k + 1; continue; } // O H3
            if (i == 3 && j == 5) { yij[k] = 0.0; k = k + 1; continue; } // O H4
            if (i == 4 && j == 5) { yij[k] = 0.0; k = k + 1; continue; } // H3 H4

            drx = x[3*i    ] - x[3*j    ];
            dry = x[3*i + 1] - x[3*j + 1];
            drz = x[3*i + 2] - x[3*j + 2];
            double dst = std::sqrt(drx*drx + dry*dry + drz*drz);

            yij[k] = std::exp(-dst / a0);
            k++;
        }
    }
}

static void make_dydr_1_2_1_2_5_intermolecular(Eigen::Ref<Eigen::MatrixXd> dydr, const double *x)
{
    double drx, dry, drz;
    size_t k = 0;

    for (size_t i = 0; i < (size_t)natoms; ++i) {
        for (size_t j = i + 1; j < (size_t)natoms; ++j) {
            drx = x[3*i    ] - x[3*j    ];
            dry = x[3*i + 1] - x[3*j + 1];
            drz = x[3*i + 2] - x[3*j + 2];

            double dst = std::sqrt(drx*drx + dry*dry + drz*drz);

            if (i == 0 && j == 1) { dydr(k, k) = 0.0; k = k + 1; continue; } // O H1
            if (i == 0 && j == 2) { dydr(k, k) = 0.0; k = k + 1; continue; } // O H2
            if (i == 1 && j == 2) { dydr(k, k) = 0.0; k = k + 1; continue; } // H1 H2
            if (i == 3 && j == 4) { dydr(k, k) = 0.0; k = k + 1; continue; } // O H3
            if (i == 3 && j == 5) { dydr(k, k) = 0.0; k = k + 1; continue; } // O H4
            if (i == 4 && j == 5) { dydr(k, k) = 0.0; k = k + 1; continue; } // H3 H4

            dydr(k, k) = -1.0/a0 * exp(-dst/a0);
            k++;
        }
    }
}

// Wrappers expected by mlp.hpp
void EVPOLY(double *y, Eigen::Ref<Eigen::RowVectorXd> p)          { evpoly_1_2_1_2_5_intermolecular(y, p); }
void EVPOLY_JAC(Eigen::Ref<Eigen::MatrixXd> jac, double *y)       { evpoly_jac_1_2_1_2_5_intermolecular(jac, y); }
void MAKE_YIJ(const double *x, double *y)                         { make_yij_1_2_1_2_5_intermolecular(x, y, natoms); }
void MAKE_DYDR(Eigen::Ref<Eigen::MatrixXd> dydr, const double *x) { make_dydr_1_2_1_2_5_intermolecular(dydr, x); }

// ---- MLP implementation (header-only library) ----

#define MLP_IMPLEMENTATION
#include "PES-IDS/mlp.hpp"

// ---- NN PES model ----

static MLPES model;

static void init_h2oh2o_nn_pes()
{
    model.init(MODEL_NPZ, natoms);
}

static double h2oh2o_nn_pes(double cart[3][6])
{
    double xyz[18] = {cart[0][0], cart[1][0], cart[2][0],
                      cart[0][1], cart[1][1], cart[2][1],
                      cart[0][2], cart[1][2], cart[2][2],
                      cart[0][3], cart[1][3], cart[2][3],
                      cart[0][4], cart[1][4], cart[2][4],
                      cart[0][5], cart[1][5], cart[2][5]};

    return model.forward(xyz);
}

static void h2oh2o_nn_pes_derivatives(double dxyz[18], double cart[3][6])
{
    double xyz[18] = {cart[0][0], cart[1][0], cart[2][0],
                      cart[0][1], cart[1][1], cart[2][1],
                      cart[0][2], cart[1][2], cart[2][2],
                      cart[0][3], cart[1][3], cart[2][3],
                      cart[0][4], cart[1][4], cart[2][4],
                      cart[0][5], cart[1][5], cart[2][5]};
    model.forward(xyz);
    memset(dxyz, 0, 18 * sizeof(double));
    model.backward(xyz, dxyz);
}

// ---- HAWAIIHYBRID interface ----

extern "C" {

void pes_init(void)
{
    init_h2oh2o_nn_pes();
}

// Minimum squared inter-monomer atom-atom distance (Bohr^2).
// Atoms closer than this trigger the repulsive wall, preventing
// NN extrapolation to spurious deep minima.
static const double RMIN_INTER_SQ = 4.0 * 4.0; // 4 Bohr ≈ 2.1 Å

static bool inter_monomer_too_close(double cart[3][6])
{
    // Monomer A: atoms 0,1,2  Monomer B: atoms 3,4,5
    for (int a = 0; a < 3; ++a) {
        for (int b = 3; b < 6; ++b) {
            double dx = cart[0][a] - cart[0][b];
            double dy = cart[1][a] - cart[1][b];
            double dz = cart[2][a] - cart[2][b];
            if (dx*dx + dy*dy + dz*dz < RMIN_INTER_SQ)
                return true;
        }
    }
    return false;
}

double pes_lab(double *q)
{
    double r_ang[9] = {q[2], q[0], q[1], q[3], q[4], q[5], q[6], q[7], q[8]};

    if (r_ang[0] <= 4.0)
        return 1e5;

    double cart[3][6];
    h2o_h2o_lab_to_cart(r_ang, cart);

    if (inter_monomer_too_close(cart))
        return 1e5;

    // h2oh2o_nn_pes returns energy in cm^-1, convert to Hartree
    return h2oh2o_nn_pes(cart) / HTOCM;
}

void dpes_lab(double *q, double *dpesdq)
{
    double r_ang[9] = {q[2], q[0], q[1], q[3], q[4], q[5], q[6], q[7], q[8]};

    if (r_ang[0] <= 4.0) {
        for (int i = 0; i < 9; ++i)
            dpesdq[i] = 0.0;
        return;
    }

    double cart[3][6];
    h2o_h2o_lab_to_cart(r_ang, cart);

    if (inter_monomer_too_close(cart)) {
        for (int i = 0; i < 9; ++i)
            dpesdq[i] = 0.0;
        return;
    }

    // Get Cartesian derivatives from NN (in cm^-1 / Bohr)
    double dxyz[18];
    h2oh2o_nn_pes_derivatives(dxyz, cart);

    Eigen::VectorXd dVdcart(18);
    for (int i = 0; i < 18; ++i)
        dVdcart(i) = dxyz[i] / HTOCM; // convert to Hartree / Bohr

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
