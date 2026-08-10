#ifndef FI_NN_PES_H
#define FI_NN_PES_H
// Same API/contract as cpp_export/fi_nn_pes.hpp -- units eV/eV-per-Bohr,
// Dalian O O H H H H atom order, Bohr, atom-major. This package has no
// fallback path: if init() succeeds, the fast generated routes are used
// unconditionally (no self-verification -- see README.md).
#ifdef __cplusplus
extern "C" {
#endif

void fi_nn_pes_init(const char* npz_path);
double fi_nn_pes_energy(const double cart18[18]);
void fi_nn_pes_energy_and_gradient(const double cart18[18], double* energy, double grad18[18]);

#ifdef __cplusplus
}
#endif
#endif  // FI_NN_PES_H
