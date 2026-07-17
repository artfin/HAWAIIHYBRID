// Sample bound H2O-Ar trajectory on the PIP-NN PES (PES-IDS/ai_pes_h2o_ar_nn_lib.cpp).
//
// Checks the analytic derivatives of the Hamiltonian against numerical ones, samples a
// bound initial condition (H < 0) with the same q_generator/p_generator/reject machinery
// the correlation drivers use, and propagates it, reporting energy conservation. No
// induced dipole surface is involved, so this exercises the PES and its derivative only.
//
// Note on the bound criterion: a pair is bound when H < 0, which presumes V -> 0 at large
// R. This PES is a raw NN extrapolation beyond its 20 a0 training range and flattens to
// -0.0826 cm-1 rather than to 0, so the true dissociation threshold sits that far below
// zero. The offset is negligible against the -140.7 cm-1 well as long as sampled energies
// are not within ~0.1 cm-1 of the threshold, which SAMPLER_RMAX below keeps them from being.

// Optionally dumps the trajectory as a multi-frame XYZ file:
//
//   ./examples/trajectory_h2o_ar.exe traj-h2o-ar.xyz
//
// Frames are written every XYZ_STRIDE steps in the frame where the H2O center of
// mass sits at the origin, so the H2O only rotates while the Ar orbits it.

#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

extern "C" {
    void pes_init(void);
    double pes_lab(double *q);
    void dpes_lab(double *q, double *dpesdq);
}

#define XYZ_STRIDE 1 

static const char *XYZ_ATOM_NAMES[4] = {"O", "H", "H", "Ar"};

// ms->intermediate_q must be up to date (call extract_q_and_write_into_ms first).
void write_xyz_frame(FILE *fp, MoleculeSystem *ms, double t, double E)
{
    double *q = ms->intermediate_q;
    double r_ang[6] = {q[2], q[0], q[1], q[3], q[4], q[5]};

    double cart[3][4];
    h2o_ar_lab_to_cart(r_ang, cart);

    fprintf(fp, "4\n");
    fprintf(fp, "t = %.1lf a.u., R = %.6lf a0, E = %.6lf cm-1\n", t, r_ang[0], E * HTOCM);

    for (size_t i = 0; i < 4; ++i) {
        fprintf(fp, "%-2s %15.8lf %15.8lf %15.8lf\n", XYZ_ATOM_NAMES[i],
                cart[0][i] * BohrToAng, cart[1][i] * BohrToAng, cart[2][i] * BohrToAng);
    }
}

const char *var_to_cstring(int n) {
    switch (n) {
        case 0:  return "Phi";
        case 1:  return "pPhi";
        case 2:  return "Theta";
        case 3:  return "pTheta";
        case 4:  return "R";
        case 5:  return "pR";
        case 6:  return "Phi1";
        case 7:  return "pPhi1";
        case 8:  return "Theta1";
        case 9:  return "pTheta1";
        case 10: return "Psi1";
        case 11: return "pPsi1";
    }

    UNREACHABLE("var_to_cstring");
}

void test_rhs(MoleculeSystem *ms, Array qp)
{
    printf("\n-----------------------------------------\n");
    printf("Testing analytic derivatives of Hamiltonian against numerical ones\n");
    printf("The derivatives are shown in the same order used in 'rhs' function.\n");

    put_qp_into_ms(ms, qp);

    N_Vector y    = make_vector(ms->QP_SIZE);
    N_Vector ydot = make_vector(ms->QP_SIZE);

    memcpy(N_VGetArrayPointer(y), qp.data, qp.n * sizeof(double));
    rhs(0.0, y, ydot, (void*) ms);

    size_t order = 6;
    Array num_derivatives = compute_numerical_rhs(ms, order);

    printf("# \t analytic \t numeric \t difference \n");
    for (size_t i = 0; i < ms->QP_SIZE; ++i) {
        printf("dot(%s): \t %.10e \t %.10e \t %.10e\n", var_to_cstring(i), NV_Ith_S(ydot, i), num_derivatives.data[i], NV_Ith_S(ydot, i) - num_derivatives.data[i]);
    }

    for (size_t i = 0; i < ms->QP_SIZE; ++i) {
        if (assert_float_is_equal_to(NV_Ith_S(ydot, i), num_derivatives.data[i], 1e-4) > 0) {
            printf("ERROR: The element (%zu) disagree!\n", i);
            exit(1);
        }
    }

    printf("-----------------------------------------\n\n\n");

    put_qp_into_ms(ms, qp);

    free_array(&num_derivatives);
    N_VDestroy(y);
    N_VDestroy(ydot);
}

int main(int argc, char *argv[])
{
    int seed = 42;

    const char *xyz_filename = (argc > 1) ? argv[1] : NULL;
    FILE *xyz_fp = NULL;

    if (xyz_filename != NULL) {
        xyz_fp = fopen(xyz_filename, "w");
        if (xyz_fp == NULL) {
            fprintf(stderr, "ERROR: could not open '%s' for writing: %s\n", xyz_filename, strerror(errno));
            exit(1);
        }
    }

    double tolerance = 1e-15; 
    size_t NSTEPS    = 1000000; 

    double MU = m_H2O * m_Ar / (m_H2O + m_Ar);
    double I1[3] = {II_H2O_A, II_H2O_B, II_H2O_C};
    MoleculeSystem ms = init_ms(MU, ROTOR, ATOM, I1, NULL, seed);

    pes_init();
    pes = pes_lab;
    dpes = dpes_lab;

    printf("H2O-Ar: reduced mass = %.10e a.u.\n", MU);

    // Global minimum of this PES, located by multistart search over the 6 lab coordinates:
    // V = -140.72 cm-1 at R = 6.93 a0. (The H2O-Ar well is known to be about -143 cm-1.)
    double Temperature = 100.0;

    CalcParams params = {};
    params.sampling_time = 200.0;
    params.ps           = PAIR_STATE_BOUND;
    params.pesmin       = -6.411647154378e-04; // Hartree = -140.72 cm-1
    params.sampler_Rmin = 4.5;                 // short-range edge of the NN training range
    params.sampler_Rmax = 12.0;                // well outside the well, but far from the -0.08 cm-1 asymptote

    // Sample a bound phase point exactly as the correlation drivers do: propose from
    // exp(-K/kT), accept with probability exp(-H/kT), then keep only H < 0.
    Array qp = create_array(ms.QP_SIZE);

    size_t attempts = 0;
    double E0;
    for (;; ++attempts) {
        q_generator(&ms, &params);
        p_generator(&ms, Temperature);

        if (reject(&ms, Temperature, params.pesmin)) continue;

        E0 = Hamiltonian(&ms);
        if (E0 < 0.0) break;
    }

    get_qp_from_ms(&ms, &qp);

    printf("Found a bound initial condition after %zu attempts\n", attempts + 1);
    printf("Initial energy: %.10e Hartree (%.4lf cm-1)\n", E0, E0 * HTOCM);
    print_array(qp);

    test_rhs(&ms, qp);

    Trajectory traj = init_trajectory(&ms, tolerance);

    E0 = Hamiltonian(&ms);
    set_initial_condition(&traj, qp);

    printf("CVODE relative tolerance: %.1e, steps: %zu\n\n", tolerance, NSTEPS);
    printf("%10s \t %12s \t %20s \t %12s\n", "t", "R", "E-E0", "V");

    double t = 0.0;
    double tout = params.sampling_time;
    double Rmin = 1e30, Rmax = -1e30, maxdE = 0.0;

    for (size_t nstep = 0; nstep < NSTEPS; ++nstep, tout += params.sampling_time)
    {
        int status = make_step(&traj, tout, &t);

        if (status > 0) {
            fprintf(stderr, "CVODE ERROR: status = %d\n", status);
            exit(1);
        }

        double E = Hamiltonian(&ms);
        extract_q_and_write_into_ms(&ms);
        double V = (*pes)(ms.intermediate_q);
        double R = ms.intermolecular_qp[IR];

        if (R < Rmin) Rmin = R;
        if (R > Rmax) Rmax = R;
        if (fabs(E - E0) > maxdE) maxdE = fabs(E - E0);

        // The trajectory is long; report a sparse sample of it.
        if (nstep % 500 == 0 || nstep == NSTEPS - 1) {
            printf("%10.1lf \t %12.10lf \t %20.15lf \t %12.6lf\n", t, R, E - E0, V * HTOCM);
        }

        if (xyz_fp != NULL && (nstep % XYZ_STRIDE == 0 || nstep == NSTEPS - 1)) {
            write_xyz_frame(xyz_fp, &ms, t, E);
        }
    }

    printf("\n--- %zu steps completed ---\n", NSTEPS);
    printf("R range visited:        %.4lf .. %.4lf a0\n", Rmin, Rmax);
    printf("max |E-E0| over run:    %.4e Hartree (%.3e cm-1)\n", maxdE, maxdE * HTOCM);
    printf("relative energy drift:  %.4e\n", maxdE / fabs(E0));

    if (xyz_fp != NULL) {
        fclose(xyz_fp);
        printf("trajectory written to:  %s (every %d steps)\n", xyz_filename, XYZ_STRIDE);
    }

    free_trajectory(&traj);
    free_ms(&ms);
    free_array(&qp);

    return 0;
}
