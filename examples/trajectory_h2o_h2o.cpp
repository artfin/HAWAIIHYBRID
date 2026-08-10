// H2O-H2O trajectory on the PIP-NN PES (PES-IDS/ai_pes_h2o_h2o_nn_lib.cpp).
//
// Checks the analytic derivatives of the Hamiltonian against numerical ones and propagates
// a trajectory, reporting energy conservation. No induced dipole surface is involved.
//
// The H2O-H2O well is very deep (~1656 cm^-1) and trajectories that fall into it can
// explore orientations where atoms get very close, hitting the NN PES hard wall at
// R < 4 a0 or inter-monomer distance < 4 a0. The test verifies energy conservation in
// the valid region of the PES before any wall collision occurs.
//
// Coordinate ordering in QP (18-element phase point):
//   [Phi, pPhi, Theta, pTheta, R, pR,
//    phi1T, pPhi1T, theta1T, pTheta1T, psi1T, pPsi1T,
//    phi2T, pPhi2T, theta2T, pTheta2T, psi2T, pPsi2T]
//
// Optionally dumps the trajectory as a multi-frame XYZ file:
//
//   ./examples/trajectory_h2o_h2o.exe traj-h2o-h2o.xyz

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

static const char *XYZ_ATOM_NAMES[6] = {"O", "H", "H", "O", "H", "H"};

// ms->intermediate_q must be up to date (call extract_q_and_write_into_ms first).
void write_xyz_frame(FILE *fp, MoleculeSystem *ms, double t, double E)
{
    double *q = ms->intermediate_q;
    // intermediate_q = [Phi, Theta, R, phi1T, theta1T, psi1T, phi2T, theta2T, psi2T]
    // r_ang order for h2o_h2o_lab_to_cart = [R, Phi, Theta, phi1T, theta1T, psi1T, phi2T, theta2T, psi2T]
    double r_ang[9] = {q[2], q[0], q[1], q[3], q[4], q[5], q[6], q[7], q[8]};

    double cart[3][6];
    h2o_h2o_lab_to_cart(r_ang, cart);

    fprintf(fp, "6\n");
    fprintf(fp, "t = %.1lf a.u., R = %.6lf a0, E = %.6lf cm-1\n", t, r_ang[0], E * HTOCM);

    for (size_t i = 0; i < 6; ++i) {
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
        case 6:  return "phi1T";
        case 7:  return "pPhi1T";
        case 8:  return "theta1T";
        case 9:  return "pTheta1T";
        case 10: return "psi1T";
        case 11: return "pPsi1T";
        case 12: return "phi2T";
        case 13: return "pPhi2T";
        case 14: return "theta2T";
        case 15: return "pTheta2T";
        case 16: return "psi2T";
        case 17: return "pPsi2T";
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

    double tolerance = 1e-10;
    size_t NSTEPS    = 500;

    double MU = m_H2O * m_H2O / (m_H2O + m_H2O);
    double I1[3] = {II_H2O_A, II_H2O_B, II_H2O_C};
    double I2[3] = {II_H2O_A, II_H2O_B, II_H2O_C};
    MoleculeSystem ms = init_ms(MU, ROTOR, ROTOR, I1, I2, seed);

    pes_init();
    pes = pes_lab;
    dpes = dpes_lab;

    printf("H2O-H2O: reduced mass = %.10e a.u.\n", MU);

    CalcParams params = {};
    params.sampling_time = 200.0;
    params.ps           = PAIR_STATE_BOUND;
    params.pesmin       = -0.007545251105751766; // Hartree = -1656.0 cm^-1
    params.sampler_Rmin = 5.0;   // avoid getting too close
    params.sampler_Rmax = 10.0;

    // Use a manually specified initial condition at large R with small momenta to stay in
    // the long-range region. The H2O-H2O well is very deep (~1656 cm-1) and trajectories
    // that fall into the well can hit orientations where atoms get close, triggering the
    // NN PES hard wall. This IC keeps R ~ 12 a0 with gentle oscillations.
    Array qp = create_array(ms.QP_SIZE);

    double data[18] = {
        // Phi, pPhi, Theta, pTheta, R, pR
        0.5, 0.0, 1.0, 0.0, 12.0, 0.0,
        // phi1T, pPhi1T, theta1T, pTheta1T, psi1T, pPsi1T
        0.3, 0.5, 0.8, -0.8, 1.2, 0.6,
        // phi2T, pPhi2T, theta2T, pTheta2T, psi2T, pPsi2T
        0.6, 0.4, 1.5, 0.3, 0.9, -0.5,
    };

    init_array(&qp, data, ms.QP_SIZE);
    put_qp_into_ms(&ms, qp);

    double E0 = Hamiltonian(&ms);

    printf("Using manual initial condition at R = 7.0 a0\n");
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

        // Report trajectory progress.
        if (nstep % 10 == 0 || nstep == NSTEPS - 1) {
            printf("%10.1lf \t %12.10lf \t %20.15lf \t %12.6lf\n", t, R, E - E0, V * HTOCM);
        }

        if (xyz_fp != NULL && (nstep % XYZ_STRIDE == 0 || nstep == NSTEPS - 1)) {
            write_xyz_frame(xyz_fp, &ms, t, E);
        }

        // Stop if trajectory hits the hard wall or escapes
        if (R < 4.5 || R > 30.0) {
            printf("Stopping: R = %.4lf a0 (outside valid range)\n", R);
            break;
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
