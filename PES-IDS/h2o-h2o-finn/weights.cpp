#include "weights.hpp"

#include <cmath>
#include <stdexcept>

#include "../cnpy.h"

namespace {

template <typename T>
std::vector<T> get_vec(cnpy::npz_t& npz, const std::string& key) {
    auto it = npz.find(key);
    if (it == npz.end()) {
        throw std::runtime_error("weights.npz missing required array '" + key + "'");
    }
    cnpy::NpyArray& arr = it->second;
    if (arr.word_size != sizeof(T)) {
        throw std::runtime_error("weights.npz array '" + key + "' has word_size=" +
                                  std::to_string(arr.word_size) + ", expected " + std::to_string(sizeof(T)));
    }
    return arr.as_vec<T>();
}

double get_scalar(cnpy::npz_t& npz, const std::string& key) { return get_vec<double>(npz, key).at(0); }
int get_scalar_i32(cnpy::npz_t& npz, const std::string& key) { return get_vec<int32_t>(npz, key).at(0); }

}  // namespace

Weights Weights::load_npz(const std::string& path) {
    cnpy::npz_t npz = cnpy::npz_load(path);

    Weights w;
    w.n_in = get_scalar_i32(npz, "n_in");
    w.n_hid = get_scalar_i32(npz, "n_hid");
    w.n_orbits = get_scalar_i32(npz, "n_orbits");

    w.W1 = get_vec<double>(npz, "W1");
    w.b1 = get_vec<double>(npz, "b1");
    w.W2 = get_vec<double>(npz, "W2");
    w.b2 = get_vec<double>(npz, "b2").at(0);

    w.in_min = get_vec<double>(npz, "in_min");
    w.in_max = get_vec<double>(npz, "in_max");
    w.out_min = get_scalar(npz, "out_min");
    w.out_max = get_scalar(npz, "out_max");

    w.radial_powers = get_vec<int32_t>(npz, "radial_powers");

    w.transform_n = get_scalar_i32(npz, "transform_n");
    const double transform_delta = get_scalar(npz, "transform_delta");
    w.transform_delta_pow_n = std::pow(transform_delta, w.transform_n);

    w.use_taper = get_scalar_i32(npz, "use_taper") != 0;
    w.taper_r_lo_bohr = get_scalar(npz, "taper_r_lo_bohr");
    w.taper_r_hi_bohr = get_scalar(npz, "taper_r_hi_bohr");

    w.bohr_to_ang = get_scalar(npz, "bohr_to_ang");

    if (static_cast<int>(w.W1.size()) != w.n_hid * w.n_in) throw std::runtime_error("W1 size mismatch");
    if (static_cast<int>(w.b1.size()) != w.n_hid) throw std::runtime_error("b1 size mismatch");
    if (static_cast<int>(w.W2.size()) != w.n_hid) throw std::runtime_error("W2 size mismatch");
    if (static_cast<int>(w.in_min.size()) != w.n_in || static_cast<int>(w.in_max.size()) != w.n_in) {
        throw std::runtime_error("in_min/in_max size mismatch");
    }
    const int expected_n_in = w.n_orbits + 5 * static_cast<int>(w.radial_powers.size());
    if (expected_n_in != w.n_in) {
        throw std::runtime_error("dimension mismatch: n_orbits + 5*radial_powers != n_in");
    }

    return w;
}
