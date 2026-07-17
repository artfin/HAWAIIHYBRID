// Shared library wrapper for the H2O-Ar PIP-NN PES
// Exposes pes_lab, dpes_lab, pes_init for HAWAIIHYBRID
//
// Model: gradient-2.npz (99 -> 64 -> 1, SiLU), trained on interaction energies and
// gradients of pes6d_Arh2o.f over R in [4.5, 20] a0 with a flexible H2O monomer.
// Permutationally invariant basis "1 2 1" (O | H H | Ar), order 4, purified;
// the generated basis/jacobian use a float ABI, so values are cast at the boundary.
//
// Input transform: all six interatomic distances enter as EXP variables,
// yij = exp(-r / a0). This must match the INTERMOLECULAR: EXP convention of the
// gradient*.yaml runs; energy.npz was trained with SWITCH-EXP6 instead and would
// need a different make_yij/make_dydr pair.
//
// Coordinate mapping:
//   HAWAIIHYBRID intermediate_q = [Phi, Theta, R, phi1T, theta1T, psi1T]
//   PES code r_ang              = [R,   Phi,   Theta, phi1T, theta1T, psi1T]
// Monomer 1 is the H2O rotor, monomer 2 is the Ar atom.

#include <cstring>

#include <Eigen/Dense>

#include "src/constants.h"
#include "src/angles_handler.hpp"

#include "PES-IDS/c_basis_1_2_1_4_purify.h"
#include "PES-IDS/c_jac_1_2_1_4_purify.h"

static const int natoms = 4;                            // O H1 H2 Ar
static const int ndist  = natoms * (natoms - 1) / 2;    // 6
static const double a0  = 2.0;                          // EXP_LAMBDA of the training runs

#define MODEL_NPZ "PES-IDS/npz/h2o-ar-pip-nn-gradient-2.npz"

// ---- PIP coordinate transforms ----

// Pair order follows combinations(4, 2): (O,H1) (O,H2) (O,Ar) (H1,H2) (H1,Ar) (H2,Ar).
static void make_yij_1_2_1_4_purify(const double *x, double *yij)
{
    size_t k = 0;
    for (int i = 0; i < natoms; ++i) {
        for (int j = i + 1; j < natoms; ++j) {
            double drx = x[3*i    ] - x[3*j    ];
            double dry = x[3*i + 1] - x[3*j + 1];
            double drz = x[3*i + 2] - x[3*j + 2];
            double dst = std::sqrt(drx*drx + dry*dry + drz*drz);

            yij[k] = std::exp(-dst / a0);
            k++;
        }
    }
}

static void make_dydr_1_2_1_4_purify(Eigen::Ref<Eigen::MatrixXd> dydr, const double *x)
{
    size_t k = 0;
    for (int i = 0; i < natoms; ++i) {
        for (int j = i + 1; j < natoms; ++j) {
            double drx = x[3*i    ] - x[3*j    ];
            double dry = x[3*i + 1] - x[3*j + 1];
            double drz = x[3*i + 2] - x[3*j + 2];
            double dst = std::sqrt(drx*drx + dry*dry + drz*drz);

            dydr(k, k) = -1.0/a0 * std::exp(-dst / a0);
            k++;
        }
    }
}

// Wrappers expected by mlp.hpp. The generated basis is float-typed, so cast on the way
// in and out; the jacobian writes only its non-zero entries and must be zeroed first.
void EVPOLY(double *y, Eigen::Ref<Eigen::RowVectorXd> p)
{
    float yf[ndist];
    for (int i = 0; i < ndist; ++i) yf[i] = static_cast<float>(y[i]);

    Eigen::RowVectorXf pf = Eigen::RowVectorXf::Zero(p.size());
    evpoly_1_2_1_4_purify(yf, pf);
    p = pf.cast<double>();
}

void EVPOLY_JAC(Eigen::Ref<Eigen::MatrixXd> jac, double *y)
{
    float yf[ndist];
    for (int i = 0; i < ndist; ++i) yf[i] = static_cast<float>(y[i]);

    Eigen::MatrixXf jacf = Eigen::MatrixXf::Zero(jac.rows(), jac.cols());
    evpoly_jac_1_2_1_4_purify(jacf, yf);
    jac = jacf.cast<double>();
}

void MAKE_YIJ(const double *x, double *y)                         { make_yij_1_2_1_4_purify(x, y); }
void MAKE_DYDR(Eigen::Ref<Eigen::MatrixXd> dydr, const double *x) { make_dydr_1_2_1_4_purify(dydr, x); }

// ---- MLP implementation (header-only library) ----

#define MLP_IMPLEMENTATION
#include "PES-IDS/mlp.hpp"

// ---- NN PES model ----

static MLPES model;

static double h2oar_nn_pes(double cart[3][4])
{
    double xyz[12] = {cart[0][0], cart[1][0], cart[2][0],
                      cart[0][1], cart[1][1], cart[2][1],
                      cart[0][2], cart[1][2], cart[2][2],
                      cart[0][3], cart[1][3], cart[2][3]};

    return model.forward(xyz);
}

static void h2oar_nn_pes_derivatives(double dxyz[12], double cart[3][4])
{
    double xyz[12] = {cart[0][0], cart[1][0], cart[2][0],
                      cart[0][1], cart[1][1], cart[2][1],
                      cart[0][2], cart[1][2], cart[2][2],
                      cart[0][3], cart[1][3], cart[2][3]};

    // `backward` reuses the activations cached by `forward`, so the two must stay paired.
    model.forward(xyz);
    memset(dxyz, 0, 12 * sizeof(double));
    model.backward(xyz, dxyz);
}

// ---- HAWAIIHYBRID interface ----

extern "C" {

void pes_init(void)
{
    model.init(MODEL_NPZ, natoms);
}

double pes_lab(double *q)
{
    double r_ang[6] = {q[2], q[0], q[1], q[3], q[4], q[5]};

    double cart[3][4];
    h2o_ar_lab_to_cart(r_ang, cart);

    // h2oar_nn_pes returns energy in cm^-1, convert to Hartree
    return h2oar_nn_pes(cart) / HTOCM;
}

void dpes_lab(double *q, double *dpesdq)
{
    double r_ang[6] = {q[2], q[0], q[1], q[3], q[4], q[5]};

    double cart[3][4];
    h2o_ar_lab_to_cart(r_ang, cart);

    // Cartesian derivatives from the NN (in cm^-1 / Bohr)
    double dxyz[12];
    h2oar_nn_pes_derivatives(dxyz, cart);

    Eigen::VectorXd dVdcart(12);
    for (int i = 0; i < 12; ++i)
        dVdcart(i) = dxyz[i] / HTOCM; // convert to Hartree / Bohr

    // Jacobian d(cart) / d(r_ang), shape 6x12
    Eigen::MatrixXd jacobi_mat = Eigen::MatrixXd::Zero(6, 12);
    h2o_ar_der_cart_by_rang(jacobi_mat, cart, r_ang);

    // dV/d(r_ang) = J * dV/d(cart)
    Eigen::Matrix<double, 6, 1> dVdrang = jacobi_mat * dVdcart;

    // Remap from r_ang order [dV/dR, dV/dPhi, dV/dTheta, ...]
    // to intermediate_q order [dV/dPhi, dV/dTheta, dV/dR, ...]
    dpesdq[0] = dVdrang(1); // dV/dPhi
    dpesdq[1] = dVdrang(2); // dV/dTheta
    dpesdq[2] = dVdrang(0); // dV/dR
    dpesdq[3] = dVdrang(3); // dV/dphi1T
    dpesdq[4] = dVdrang(4); // dV/dtheta1T
    dpesdq[5] = dVdrang(5); // dV/dpsi1T
}

} // extern "C"
