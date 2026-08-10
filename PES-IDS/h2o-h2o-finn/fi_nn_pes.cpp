#include "fi_nn_pes.hpp"

#include <cstdio>
#include <cstdlib>

#include "fast_pes_pipeline.hpp"
#include "weights.hpp"

namespace {

Weights g_weights;
bool g_initialized = false;

void require_initialized() {
    if (!g_initialized) {
        std::fprintf(stderr, "fi_nn_pes: fi_nn_pes_init() must be called before any energy call\n");
        std::abort();
    }
}

}  // namespace

extern "C" void fi_nn_pes_init(const char* npz_path) {
    g_weights = Weights::load_npz(npz_path);

    // No runtime AD/generic fallback exists in this package (see README.md)
    // -- fail loudly at init if the generated files don't match this
    // checkpoint's dimensions, rather than silently compute wrong answers.
    if (g_weights.n_orbits != fi_nn_pes::N_ORBITS_GENERATED) {
        std::fprintf(stderr,
                      "fi_nn_pes (fast): n_orbits=%d != N_ORBITS_GENERATED=%d -- "
                      "orbit_basis_generated.hpp is stale relative to this checkpoint. "
                      "Regenerate it in cpp_export/, re-verify there, then re-sync into cpp_export_fast/.\n",
                      g_weights.n_orbits, fi_nn_pes::N_ORBITS_GENERATED);
        std::abort();
    }
    if (g_weights.n_in != fi_nn_pes::N_FEATURES_GENERATED) {
        std::fprintf(stderr,
                      "fi_nn_pes (fast): n_in=%d != N_FEATURES_GENERATED=%d -- "
                      "features_jacobian_generated.hpp is stale relative to this checkpoint. "
                      "Regenerate it in cpp_export/, re-verify there, then re-sync into cpp_export_fast/.\n",
                      g_weights.n_in, fi_nn_pes::N_FEATURES_GENERATED);
        std::abort();
    }

    g_initialized = true;
}

extern "C" double fi_nn_pes_energy(const double cart18[18]) {
    require_initialized();
    return fi_nn_pes::pip_pes_energy_fast(cart18, g_weights);
}

extern "C" void fi_nn_pes_energy_and_gradient(const double cart18[18], double* energy, double grad18[18]) {
    require_initialized();
    *energy = fi_nn_pes::pip_pes_energy_and_gradient_fast(cart18, g_weights, grad18);
}
