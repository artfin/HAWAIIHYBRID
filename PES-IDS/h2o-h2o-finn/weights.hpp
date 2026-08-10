#pragma once
// Trimmed checkpoint container for the FAST-ONLY package (cpp_export_fast/).
// Drops orbit_exps/orbit_term_counts/max_exponent -- those only feed
// evaluate_orbit_basis_generic, which this package deliberately excludes
// (see README.md: this package is the two already-measured-fastest routes
// only, no generic/AD fallback). Everything else matches cpp_export/
// weights.hpp exactly, loaded from the SAME weights.npz.
#include <cstdint>
#include <string>
#include <vector>

constexpr int MAX_HIDDEN = 128;
constexpr int MAX_RADIAL_POWERS = 8;
constexpr int MAX_FEATURES = 512 + 5 * MAX_RADIAL_POWERS;

struct Weights {
    int n_in = 0;
    int n_hid = 0;
    int n_orbits = 0;

    std::vector<double> W1;  // (n_hid, n_in) row-major
    std::vector<double> b1;  // (n_hid)
    std::vector<double> W2;  // (n_hid)
    double b2 = 0.0;

    std::vector<double> in_min;
    std::vector<double> in_max;
    double out_min = 0.0;
    double out_max = 0.0;

    std::vector<int32_t> radial_powers;

    int transform_n = 2;
    double transform_delta_pow_n = 4.0;

    bool use_taper = true;
    double taper_r_lo_bohr = 25.0;
    double taper_r_hi_bohr = 35.0;

    double bohr_to_ang = 0.529177212;

    static Weights load_npz(const std::string& path);
};
