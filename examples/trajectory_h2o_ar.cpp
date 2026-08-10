// Sample bound H2O-Ar trajectory on the PIP-NN PES (PES-IDS/ai_pes_h2o_ar_nn_lib.cpp).
//
// Checks the analytic derivatives of the Hamiltonian against numerical ones, samples a
// bound initial condition (H < 0) with the same q_generator/p_generator/reject machinery
// the correlation drivers use, and propagates it, reporting energy conservation. No
// induced dipole surface is involved, so this exercises the PES and its derivative only.
//
// The PES wrapper removes the NN's constant large-R residual, so H < 0 uses the physical
// dissociation zero directly.

// Optionally dumps the trajectory as a multi-frame XYZ file:
//
//   ./examples/trajectory_h2o_ar.exe traj.xyz diagnostics.csv [steps] [stride]
//       [bound|nonnegative|all] [seed] [trajectories] [energy_max_cm-1]
//       [temperature_K] [energy_min_cm-1] [restart_file|-] [time_offset_au]
//       [absorbing_escape_radius_bohr]
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

typedef struct {
    const char *trajectory_id;
    FILE *stream;
    size_t stride;
    size_t step;
    bool has_outer_barrier;
    double G_hartree;
    double R_barrier_bohr;
    double temperature_K;
    double time_offset_au;
    double last_written_time_au;
} TrajectoryMetadata;

static double vector_norm(const double v[3])
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static void write_diagnostic_header(FILE *stream, const MoleculeSystem *ms)
{
    fprintf(stream, "trajectory_id,t_au,temperature_K,R_bohr,E_total_hartree,V_hartree,K_hartree,pR_au,");
    fprintf(stream, "Lx_au,Ly_au,Lz_au,L_norm_au,Jrot_x_au,Jrot_y_au,Jrot_z_au,Jrot_norm_au,");
    fprintf(stream, "Jtotal_x_au,Jtotal_y_au,Jtotal_z_au,Jtotal_norm_au,");
    fprintf(stream, "has_outer_barrier,G_hartree,R_barrier_bohr");
    for (size_t i = 0; i < ms->QP_SIZE; ++i) fprintf(stream, ",qp_%zu", i);
    fputc('\n', stream);
}

static void write_diagnostic_row(const Trajectory *traj, double t, TrajectoryMetadata *metadata)
{
    t += metadata->time_offset_au;
    if (t == metadata->last_written_time_au) return;
    metadata->last_written_time_au = t;
    MoleculeSystem *ms = traj->ms;
    extract_q_and_write_into_ms(ms);
    double V = pes(ms->intermediate_q), K = kinetic_energy(ms);
    double L[3], Jrot[3], Jtotal[3];
    j_orbital(ms, L);
    j_monomer(&ms->m1, Jrot);
    for (size_t k = 0; k < 3; ++k) Jtotal[k] = L[k] + Jrot[k];

    FILE *stream = metadata->stream;
    fprintf(stream, "%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,",
            metadata->trajectory_id, t, metadata->temperature_K,
            ms->intermolecular_qp[IR], V + K, V, K,
            ms->intermolecular_qp[IPR]);
    fprintf(stream, "%.17g,%.17g,%.17g,%.17g,", L[0], L[1], L[2], vector_norm(L));
    fprintf(stream, "%.17g,%.17g,%.17g,%.17g,", Jrot[0], Jrot[1], Jrot[2], vector_norm(Jrot));
    fprintf(stream, "%.17g,%.17g,%.17g,%.17g,", Jtotal[0], Jtotal[1], Jtotal[2], vector_norm(Jtotal));
    fprintf(stream, "%d,%.17g,%.17g", metadata->has_outer_barrier ? 1 : 0,
            metadata->G_hartree, metadata->R_barrier_bohr);
    Array qp = create_array(ms->QP_SIZE);
    get_qp_from_ms(ms, &qp);
    for (size_t i = 0; i < qp.n; ++i) fprintf(stream, ",%.17g", qp.data[i]);
    fputc('\n', stream);
    free_array(&qp);
}

static void write_diagnostic_sample(const Trajectory *traj, double t, void *user_data)
{
    TrajectoryMetadata *metadata = (TrajectoryMetadata *) user_data;
    if (metadata->step++ % metadata->stride != 0) return;
    write_diagnostic_row(traj, t, metadata);
}

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
    int seed = (argc > 6) ? atoi(argv[6]) : 42;
    const char *energy_class = (argc > 5) ? argv[5] : "bound";
    if (strcmp(energy_class, "bound") != 0 && strcmp(energy_class, "nonnegative") != 0 &&
        strcmp(energy_class, "all") != 0) {
        fprintf(stderr, "ERROR: energy class must be bound, nonnegative, or all\n");
        exit(1);
    }

    const char *xyz_filename = (argc > 1 && strcmp(argv[1], "-") != 0) ? argv[1] : NULL;
    const char *csv_filename = (argc > 2) ? argv[2] : NULL;
    FILE *xyz_fp = NULL;
    FILE *csv_fp = NULL;

    if (xyz_filename != NULL) {
        xyz_fp = fopen(xyz_filename, "w");
        if (xyz_fp == NULL) {
            fprintf(stderr, "ERROR: could not open '%s' for writing: %s\n", xyz_filename, strerror(errno));
            exit(1);
        }
    }
    if (csv_filename != NULL) {
        csv_fp = fopen(csv_filename, "w");
        if (csv_fp == NULL) {
            fprintf(stderr, "ERROR: could not open '%s' for writing: %s\n", csv_filename, strerror(errno));
            exit(1);
        }
    }

    double tolerance = 1e-15; 
    size_t NSTEPS = (argc > 3) ? strtoull(argv[3], NULL, 10) : 1000000;
    size_t output_stride = (argc > 4) ? strtoull(argv[4], NULL, 10) : XYZ_STRIDE;
    size_t ntrajectories = (argc > 7) ? strtoull(argv[7], NULL, 10) : 1;
    double energy_max = (argc > 8) ? atof(argv[8]) / HTOCM : INFINITY;
    double Temperature = (argc > 9) ? atof(argv[9]) : 296.0;
    double energy_min = (argc > 10) ? atof(argv[10]) / HTOCM : -INFINITY;
    const char *restart_filename = (argc > 11 && strcmp(argv[11], "-") != 0) ? argv[11] : NULL;
    double time_offset = (argc > 12) ? atof(argv[12]) : 0.0;
    double escape_radius = (argc > 13) ? atof(argv[13]) : INFINITY;
    if (NSTEPS == 0 || output_stride == 0 || ntrajectories == 0) {
        fprintf(stderr, "ERROR: steps, stride, and trajectories must be positive\n");
        exit(1);
    }

    double MU = m_H2O * m_Ar / (m_H2O + m_Ar);
    double I1[3] = {II_H2O_A, II_H2O_B, II_H2O_C};
    MoleculeSystem ms = init_ms(MU, ROTOR, ATOM, I1, NULL, seed);

    pes_init();
    pes = pes_lab;
    dpes = dpes_lab;

    printf("H2O-Ar: reduced mass = %.10e a.u.\n", MU);

    // Global minimum of this PES, located by multistart search over the 6 lab coordinates:
    // V = -140.72 cm-1 at R = 6.93 a0. (The H2O-Ar well is known to be about -143 cm-1.)
    CalcParams params = {};
    params.sampling_time = 200.0;
    params.ps           = PAIR_STATE_BOUND;
    params.pesmin       = -6.411647154378e-04; // Hartree = -140.72 cm-1
    params.sampler_Rmin = 4.5;                 // short-range edge of the NN training range
    params.sampler_Rmax = 12.0;                // well outside the well, but far from the -0.08 cm-1 asymptote

    Array qp = create_array(ms.QP_SIZE);
    if (csv_fp != NULL) write_diagnostic_header(csv_fp, &ms);
    FILE *restart_fp = NULL;
    if (restart_filename != NULL) {
        restart_fp = fopen(restart_filename, "r");
        if (restart_fp == NULL) {
            fprintf(stderr, "ERROR: could not open restart file '%s': %s\n", restart_filename, strerror(errno));
            exit(1);
        }
    }

    size_t total_attempts = 0;
    for (size_t trajectory_index = 0; trajectory_index < ntrajectories; ++trajectory_index) {
        size_t attempts = 0;
        double E0;
        char trajectory_id[128];
        if (restart_fp != NULL) {
            if (fscanf(restart_fp, "%127s", trajectory_id) != 1) {
                fprintf(stderr, "ERROR: restart file ended before trajectory %zu\n", trajectory_index);
                exit(1);
            }
            for (size_t i = 0; i < qp.n; ++i) {
                if (fscanf(restart_fp, "%lf", &qp.data[i]) != 1) {
                    fprintf(stderr, "ERROR: malformed restart phase point for '%s'\n", trajectory_id);
                    exit(1);
                }
            }
            put_qp_into_ms(&ms, qp);
            E0 = Hamiltonian(&ms);
        } else {
            for (;; ++attempts) {
                q_generator(&ms, &params);
                p_generator(&ms, Temperature);
                if (reject(&ms, Temperature, params.pesmin)) continue;
                E0 = Hamiltonian(&ms);
                if (strcmp(energy_class, "bound") == 0 && E0 >= 0.0) continue;
                if (strcmp(energy_class, "nonnegative") == 0 && E0 < 0.0) continue;
                if (E0 < energy_min) continue;
                if (E0 > energy_max) continue;
                break;
            }
            get_qp_from_ms(&ms, &qp);
            snprintf(trajectory_id, sizeof(trajectory_id), "h2o_ar_%s_seed_%d_%06zu",
                     energy_class, seed, trajectory_index);
        }
        total_attempts += attempts + 1;
        if (trajectory_index == 0) test_rhs(&ms, qp);

        Trajectory traj = init_trajectory(&ms, tolerance);
        set_initial_condition(&traj, qp);
        TrajectoryMetadata metadata = {
            .trajectory_id = trajectory_id, .stream = csv_fp,
            .stride = output_stride, .step = 0,
            .has_outer_barrier = false, .G_hartree = NAN, .R_barrier_bohr = NAN,
            .temperature_K = Temperature,
            .time_offset_au = time_offset, .last_written_time_au = -INFINITY,
        };
        if (csv_fp != NULL) {
            write_diagnostic_sample(&traj, 0.0, &metadata);
            trajectory_set_step_observer(&traj, write_diagnostic_sample, &metadata);
        }

        double t = 0.0, tout = params.sampling_time;
        double Rmin = ms.intermolecular_qp[IR], Rmax = Rmin, maxdE = 0.0;
        for (size_t nstep = 0; nstep < NSTEPS; ++nstep, tout += params.sampling_time) {
            int status = make_step(&traj, tout, &t);
            if (status > 0) {
                fprintf(stderr, "CVODE ERROR: status = %d, trajectory = %s\n", status, trajectory_id);
                exit(1);
            }
            double E = Hamiltonian(&ms), R = ms.intermolecular_qp[IR];
            Rmin = fmin(Rmin, R); Rmax = fmax(Rmax, R); maxdE = fmax(maxdE, fabs(E-E0));
            if (xyz_fp != NULL && (nstep % output_stride == 0 || nstep == NSTEPS-1))
                write_xyz_frame(xyz_fp, &ms, t, E);
            if (R >= escape_radius && ms.intermolecular_qp[IPR] > 0.0) {
                if (csv_fp != NULL) write_diagnostic_row(&traj, t, &metadata);
                break;
            }
        }
        printf("[%zu/%zu] %s E=%.4f cm-1 R=%.3f..%.3f max_dE=%.3e cm-1 attempts=%zu\n",
               trajectory_index+1, ntrajectories, trajectory_id, E0*HTOCM, Rmin, Rmax,
               maxdE*HTOCM, attempts+1);
        free_trajectory(&traj);
    }
    printf("Completed %zu trajectories after %zu total proposals\n", ntrajectories, total_attempts);
    if (restart_fp != NULL) fclose(restart_fp);

    if (xyz_fp != NULL) {
        fclose(xyz_fp);
        printf("trajectory written to:  %s (every %zu steps)\n", xyz_filename, output_stride);
    }
    if (csv_fp != NULL) {
        fclose(csv_fp);
        printf("trajectory diagnostics written to: %s (every %zu steps)\n", csv_filename, output_stride);
    }

    free_ms(&ms);
    free_array(&qp);

    return 0;
}
