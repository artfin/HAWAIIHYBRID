#define USE_MPI
#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

#include "ai_pes_co2ar.h"

#include "ai_ids_co2ar.hpp"
static AI_IDS_co2_ar co2_ar_ids;

double pes(double *q) {
    static double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);
    return pes_co2ar(qmol[0], qmol[4]);
} 

void dpes(double *q, double *dpesdq) {
    static Eigen::Matrix<double, 5, 5> jac;
    static Eigen::Matrix<double, 5, 1> derivatives_mol, derivatives_lab; 
    static double qmol[5];
    
    jac.setZero();
    linear_molecule_atom_lab_to_mol(q, qmol);
    linear_molecule_atom_Jacobi_mol_by_lab(jac, q, qmol);  
    
    double dR, dTheta;
    dpes_co2ar(qmol[0], qmol[4], &dR, &dTheta); // [R, THETAM] -> [dpes_dR, dpes_dTheta]

    // [dpes_dR, 0, 0, 0, dpes_dTheta]
    derivatives_mol(0) = dR; 
    derivatives_mol(4) = dTheta; 
    
    derivatives_lab = jac * derivatives_mol;
    Eigen::VectorXd::Map(dpesdq, 5) = derivatives_lab;
}

void dipole_lab(double *q, double diplab[3]) {
    double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);
    
    double dipmol[3];
    co2_ar_ids.dipole_vector(qmol[0], qmol[4], dipmol); 
    
    double sinphiem, cosphiem;
    double sinthetaem, costhetaem;
    double sinpsiem, cospsiem;

    sincos(qmol[1], &sinphiem, &cosphiem);
    sincos(qmol[2], &sinthetaem, &costhetaem);
    sincos(qmol[3], &sinpsiem, &cospsiem);

    Sz_filler(Sphiem, sinphiem, cosphiem);
    Sx_filler(Sthetaem, sinthetaem, costhetaem);
    Sz_filler(Spsiem, sinpsiem, cospsiem);
       
    Eigen::Vector3d dipmol_eig = Eigen::Map<Eigen::Vector3d>(dipmol, 3);
    Eigen::Vector3d diplab_eig = Sphiem.transpose() * Sthetaem.transpose() * Spsiem.transpose() * dipmol_eig; 
   
    diplab[0] = diplab_eig(0); 
    diplab[1] = diplab_eig(1); 
    diplab[2] = diplab_eig(2); 
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    MPI_Context ctx = {};
    MPI_Comm_size(MPI_COMM_WORLD, &ctx.size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &ctx.rank); 

    int seed = 42;

    double MU = m_CO2 * m_Ar / (m_CO2 + m_Ar); 
    double I1[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL, seed);

    init_pes();
    co2_ar_ids.init();
    dipole = dipole_lab;

    double tolerance = 1e-12;
    Trajectory traj = init_trajectory(&ms, tolerance);
   
    CalcParams params = {};
    params.sampler_Rmin                     = 4.0;
    params.sampler_Rmax                     = 40.0;
    params.total_trajectories               = 10;
    params.sampling_time                    = 200.0;
    params.MaxTrajectoryLength              = 65536;
    params.Rcut                             = 40.0;
    params.partial_partition_function_ratio = 1.0;

    double Temperature = 300.0;
    
    CFnc cf = calculate_correlation(ctx, &ms, &params, Temperature);

    if (ctx.rank == 0) {
        printf("\n\n");
        printf("Correlation function is calculated. The initial values are:\n");
        for (size_t i = 0; i < 10; ++i) {
            printf("%.2f   %.10e\n", cf.t[i], cf.vals[i]);
        }
    }

    free_trajectory(&traj);
    free_ms(&ms);

    MPI_Finalize();

    return 0;
}


