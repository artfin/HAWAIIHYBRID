#pragma once
// The fastest-route-only PES pipeline. Contains ONLY the two paths measured
// as fastest in cpp_export/ (lab_notes.md sec.37.1-37.4): the code-generated
// value-only orbit basis for energy queries, and the code-generated full
// analytic Jacobian for energy+gradient queries -- both chained with the
// cache-locality-fixed MLP forward/backward. No AD (`Var`/`Tape`), no
// generic runtime-loop orbit-basis fallback, no self-verification: this
// package TRUSTS that orbit_basis_generated.hpp/features_jacobian_generated.hpp
// were already validated upstream in cpp_export/ (which has the full AD/
// generic ground-truth machinery and verify_against_pytorch.py). See
// README.md for the sync workflow and why this tradeoff is safe in practice.
//
// Atom order: native Dalian O O H H H H, Bohr, atom-major xyz -- identical
// contract to cpp_export/fi_nn_pes.hpp.
#include <cmath>

#include "features_jacobian_generated.hpp"
#include "orbit_basis_generated.hpp"
#include "weights.hpp"

namespace fi_nn_pes {

constexpr int PIP_ORDER[6] = {0, 2, 3, 1, 4, 5};
constexpr int INTER_PAIRS[9][2] = {{0, 3}, {0, 4}, {0, 5}, {1, 3}, {1, 4}, {1, 5}, {2, 3}, {2, 4}, {2, 5}};
constexpr double PI = 3.14159265358979323846;

inline double dot3(const double a[3], const double b[3]) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline void sub3(const double a[3], const double b[3], double out[3]) {
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
}
inline void normalize3(double v[3]) {
    const double inv_norm = 1.0 / std::sqrt(dot3(v, v));
    v[0] *= inv_norm;
    v[1] *= inv_norm;
    v[2] *= inv_norm;
}

// 4-way manually-unrolled dot product -- lets the compiler auto-vectorize
// without needing -ffast-math (explicit reassociation done once in source,
// not compiler-driven). See cpp_export/README.md's "MLP backward pass had a
// transposed-matrix cache-locality bug" section for the measurement behind
// this and the loop-order choice below.
inline double dot_unrolled(const double* a, const double* b, int n) {
    double acc0 = 0.0, acc1 = 0.0, acc2 = 0.0, acc3 = 0.0;
    int j = 0;
    for (; j + 4 <= n; j += 4) {
        acc0 += a[j] * b[j];
        acc1 += a[j + 1] * b[j + 1];
        acc2 += a[j + 2] * b[j + 2];
        acc3 += a[j + 3] * b[j + 3];
    }
    double acc = (acc0 + acc1) + (acc2 + acc3);
    for (; j < n; ++j) acc += a[j] * b[j];
    return acc;
}

inline double mlp_forward(const double X[], const Weights& w) {
    double Xs[MAX_FEATURES];
    for (int j = 0; j < w.n_in; ++j) Xs[j] = 2.0 * (X[j] - w.in_min[j]) / (w.in_max[j] - w.in_min[j]) - 1.0;
    double hidden[MAX_HIDDEN];
    for (int h = 0; h < w.n_hid; ++h) {
        const double* row = &w.W1[static_cast<std::size_t>(h) * w.n_in];
        const double acc = w.b1[h] + dot_unrolled(row, Xs, w.n_in);
        hidden[h] = acc / std::sqrt(1.0 + acc * acc);
    }
    double out = w.b2;
    for (int h = 0; h < w.n_hid; ++h) out += w.W2[h] * hidden[h];
    return (out + 1.0) * (w.out_max - w.out_min) * 0.5 + w.out_min;
}

inline double mlp_forward_backward(const double X[], const Weights& w, double dE_dX_out[MAX_FEATURES]) {
    double Xs[MAX_FEATURES];
    for (int j = 0; j < w.n_in; ++j) Xs[j] = 2.0 * (X[j] - w.in_min[j]) / (w.in_max[j] - w.in_min[j]) - 1.0;

    double pre[MAX_HIDDEN], hidden[MAX_HIDDEN];
    for (int h = 0; h < w.n_hid; ++h) {
        const double* row = &w.W1[static_cast<std::size_t>(h) * w.n_in];
        const double acc = w.b1[h] + dot_unrolled(row, Xs, w.n_in);
        pre[h] = acc;
        hidden[h] = acc / std::sqrt(1.0 + acc * acc);
    }

    double out = w.b2;
    for (int h = 0; h < w.n_hid; ++h) out += w.W2[h] * hidden[h];
    const double nn_val = (out + 1.0) * (w.out_max - w.out_min) * 0.5 + w.out_min;

    const double d_out = (w.out_max - w.out_min) * 0.5;
    double d_pre[MAX_HIDDEN];
    for (int h = 0; h < w.n_hid; ++h) {
        const double v = pre[h];
        const double s = std::sqrt(1.0 + v * v);
        d_pre[h] = (w.W2[h] * d_out) / (s * s * s);
    }

    // h-outer/j-inner: contiguous W1 access (see module comment).
    for (int j = 0; j < w.n_in; ++j) dE_dX_out[j] = 0.0;
    for (int h = 0; h < w.n_hid; ++h) {
        const double dph = d_pre[h];
        const double* row = &w.W1[static_cast<std::size_t>(h) * w.n_in];
        for (int j = 0; j < w.n_in; ++j) dE_dX_out[j] += dph * row[j];
    }
    for (int j = 0; j < w.n_in; ++j) dE_dX_out[j] *= (2.0 / (w.in_max[j] - w.in_min[j]));

    return nn_val;
}

// --- energy-only: fastest measured route (generated orbit values + generic
// angle features, which are cheap enough not to be the bottleneck) ---
inline double pip_pes_energy_fast(const double cart[18], const Weights& w) {
    double c[6][3];
    for (int k = 0; k < 6; ++k) {
        const int a = PIP_ORDER[k];
        for (int d = 0; d < 3; ++d) c[k][d] = cart[3 * a + d] * w.bohr_to_ang;
    }

    double y9[9];
    for (int p = 0; p < 9; ++p) {
        const int a = INTER_PAIRS[p][0], b = INTER_PAIRS[p][1];
        double d[3];
        sub3(c[a], c[b], d);
        const double r = std::sqrt(dot3(d, d));
        double rn = r;
        for (int k = 1; k < w.transform_n; ++k) rn *= r;  // r^transform_n, never pow(r,0) (n>=1 always)
        y9[p] = 1.0 / (rn + w.transform_delta_pow_n);
    }
    const double Roo_T = y9[0];

    double X[MAX_FEATURES];
    evaluate_orbit_basis_generated(y9, X);

    double uR[3];
    sub3(c[3], c[0], uR);
    normalize3(uR);
    double uA1[3];
    sub3(c[1], c[0], uA1);
    normalize3(uA1);
    double uA2[3];
    sub3(c[2], c[0], uA2);
    normalize3(uA2);
    double uB1[3];
    sub3(c[4], c[3], uB1);
    normalize3(uB1);
    double uB2[3];
    sub3(c[5], c[3], uB2);
    normalize3(uB2);
    double negUR[3] = {-uR[0], -uR[1], -uR[2]};
    const double cosA1 = dot3(uA1, uR), cosA2 = dot3(uA2, uR);
    const double cosB1 = dot3(uB1, negUR), cosB2 = dot3(uB2, negUR);
    const double cross = dot3(uA1, uB1) + dot3(uA1, uB2) + dot3(uA2, uB1) + dot3(uA2, uB2);
    const double ang[5] = {cosA1 + cosA2, cosA1 * cosA1 + cosA2 * cosA2, cosB1 + cosB2,
                            cosB1 * cosB1 + cosB2 * cosB2, cross};

    int col = N_ORBITS_GENERATED;
    for (std::size_t pi = 0; pi < w.radial_powers.size(); ++pi) {
        const int p = w.radial_powers[pi];
        double rp = 1.0;
        for (int k = 0; k < p; ++k) rp *= Roo_T;
        for (int j = 0; j < 5; ++j) X[col++] = ang[j] * rp;
    }

    double nn_val = mlp_forward(X, w);

    if (w.use_taper) {
        const double dRx = cart[0] - cart[3], dRy = cart[1] - cart[4], dRz = cart[2] - cart[5];
        const double R = std::sqrt(dRx * dRx + dRy * dRy + dRz * dRz);
        double t = (R - w.taper_r_lo_bohr) / (w.taper_r_hi_bohr - w.taper_r_lo_bohr);
        t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
        const double taper = 1.0 - 0.5 * (1.0 - std::cos(PI * t));
        nn_val *= taper;
    }
    return nn_val;
}

// --- energy+gradient: fastest measured route (fully analytic Jacobian) ---
inline double pip_pes_energy_and_gradient_fast(const double cart[18], const Weights& w, double grad18[18]) {
    double X[N_FEATURES_GENERATED];
    double J[N_FEATURES_GENERATED][18];
    evaluate_features_and_jacobian_generated(cart, X, J);

    double dE_dX[MAX_FEATURES];
    const double nn_val = mlp_forward_backward(X, w, dE_dX);

    double grad_nn[18] = {0.0};
    for (int k = 0; k < w.n_in; ++k) {
        const double d = dE_dX[k];
        if (d == 0.0) continue;
        for (int m = 0; m < 18; ++m) grad_nn[m] += d * J[k][m];
    }

    if (!w.use_taper) {
        for (int m = 0; m < 18; ++m) grad18[m] = grad_nn[m];
        return nn_val;
    }

    const double dRx = cart[0] - cart[3], dRy = cart[1] - cart[4], dRz = cart[2] - cart[5];
    const double R = std::sqrt(dRx * dRx + dRy * dRy + dRz * dRz);
    const double t_raw = (R - w.taper_r_lo_bohr) / (w.taper_r_hi_bohr - w.taper_r_lo_bohr);

    double taper, dtaper_dR;
    if (t_raw <= 0.0) {
        taper = 1.0;
        dtaper_dR = 0.0;
    } else if (t_raw >= 1.0) {
        taper = 0.0;
        dtaper_dR = 0.0;
    } else {
        taper = 1.0 - 0.5 * (1.0 - std::cos(PI * t_raw));
        const double dt_dR = 1.0 / (w.taper_r_hi_bohr - w.taper_r_lo_bohr);
        dtaper_dR = -0.5 * std::sin(PI * t_raw) * PI * dt_dR;
    }

    const double energy = taper * nn_val;
    const double inv_R = (R > 0.0) ? (1.0 / R) : 0.0;
    const double dR_dcart[6] = {dRx * inv_R,  dRy * inv_R,  dRz * inv_R,
                                 -dRx * inv_R, -dRy * inv_R, -dRz * inv_R};
    for (int d = 0; d < 6; ++d) grad18[d] = taper * grad_nn[d] + nn_val * dtaper_dR * dR_dcart[d];
    for (int m = 6; m < 18; ++m) grad18[m] = taper * grad_nn[m];
    return energy;
}

}  // namespace fi_nn_pes
