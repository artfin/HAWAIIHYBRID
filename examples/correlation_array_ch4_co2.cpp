#define USE_MPI
#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

#include "ai_pes_ch4_co2.h"
#include "ai_ids_ch4_co2.hpp"

#define NLAB 8
#define NKAL 5

double pes(double *qlab) {
    static double qkal[NKAL];
    CH4_linear_molecule_lab_to_kal(qlab, qkal); 
    return pes_ch4co2(qkal[0], qkal[1], qkal[2], qkal[3], qkal[4]); 
} 

void dpes(double *qlab, double *dpesdq) {
    static Eigen::Matrix<double, NLAB, NKAL> jac;
    static Eigen::Matrix<double, NKAL, 1> derivatives_kal;
    static Eigen::Matrix<double, NLAB, 1> derivatives_lab; 
    static double qkal[NKAL];

    jac.setZero();
    CH4_linear_molecule_lab_to_kal(qlab, qkal);
    CH4_linear_molecule_Jacobi_kal_by_lab(jac, qlab, qkal);  

    double dR, dphi1, dphi2, dtheta1, dtheta2; 
    dpes_ch4co2(qkal[0], qkal[1], qkal[2], qkal[3], qkal[4], &dR, &dphi1, &dtheta1, &dphi2, &dtheta2); 

    derivatives_kal(0) = dR; 
    derivatives_kal(1) = dphi1;
    derivatives_kal(2) = dtheta1;
    derivatives_kal(3) = dphi2;
    derivatives_kal(4) = dtheta2;

    derivatives_lab = jac * derivatives_kal;

    Eigen::VectorXd::Map(dpesdq, NLAB) = derivatives_lab;
}

void dipole_lab(double *qlab, double diplab[3]) {
    double qkal[NKAL];
    CH4_linear_molecule_lab_to_kal(qlab, qkal);
    
    double dipkal[3];
    dipole_vector(qkal, dipkal);

    double sinphi, cosphi;
    double sintheta, costheta;
    double sinpsi, cospsi;

    sincos(qlab[3], &sinphi, &cosphi);
    sincos(qlab[4], &sintheta, &costheta);
    sincos(qlab[5], &sinpsi, &cospsi);

    Sz_filler(Sphi1t, sinphi, cosphi);
    Sx_filler(Stheta1t, sintheta, costheta);
    Sz_filler(Spsi1t, sinpsi, cospsi);
      
    Eigen::Vector3d dipkal_eig = Eigen::Map<Eigen::Vector3d>(dipkal, 3);
    Eigen::Vector3d diplab_eig = Sphi1t.transpose() * Stheta1t.transpose() * Spsi1t.transpose() * dipkal_eig; 
   
    diplab[0] = diplab_eig(0); 
    diplab[1] = diplab_eig(1); 
    diplab[2] = diplab_eig(2); 
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    INIT_WRANK;
    INIT_WSIZE;

    int seeds[] = {50377445, 37183603, 13442809, 79487979, 22867865, 65728383, 35191887, 86692103, 71931368, 49143611, 47120112, 20335931, 94038811, 67455850, 84879044, 64511517, 23926327, 18913047, 38811051, 40524201, 23443742, 69569457, 73376892, 37590224, 63813629, 24265718, 72679138, 67977324, 10974897, 87328344, 21116964, 30240032, 77558402, 77890020, 99138094, 44532637, 43714282, 96995088, 38135587, 18821238, 72311571, 58108907, 15028065, 93134752, 90678150, 24939378, 13871796, 11189243, 12183134, 42314286, 58177735, 22031114, 38428528, 54438739, 31300835, 11735683, 64433370, 60589114, 99641636, 76404246, 63249860, 53885845, 11459935, 10866339, 28665977, 85378357, 20730314, 54747850, 35369175, 80381443, 75108215, 27845484, 88480851, 88869525, 97779085, 87433123, 47820043, 30924110, 19541830, 13611664, 91726323, 40436189, 76529728, 61510199, 55802727, 51413505, 80876909, 56875599, 33484555, 95684681, 82409850, 17629235, 21195855, 31645930, 35336321, 24600932, 92287561, 98020917, 24661211, 71793140};
    assert(_wsize <= (int) (sizeof(seeds)/sizeof(int)) && "ASSERT: more seeds should be added\n");

    double MU = m_CH4 * m_CO2 / (m_CH4 + m_CO2); 
    double I1[3] = {II_CH4, II_CH4, II_CH4};
    double I2[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, ROTOR, LINEAR_MOLECULE, I1, I2, seeds[_wrank]);

    init_pes();
    init_ids();
    dipole = dipole_lab;

    double tolerance = 1e-12;
    Trajectory traj = init_trajectory(&ms, tolerance);
    
    CalcParams params = {};
    params.ps                               = FREE_AND_METASTABLE;
    params.sampler_Rmin                     = 4.751;
    params.sampler_Rmax                     = 40.0;
    params.niterations                      = 100;
    params.total_trajectories               = 10000000; 
    params.cvode_tolerance                  = 1e-12;
    params.sampling_time                    = 200.0;
    params.MaxTrajectoryLength              = 65536;
    params.Rcut                             = 40.0;
    params.initialM0_npoints                = 20000000;
    params.initialM2_npoints                = 20000000;
    params.pesmin                           = -342.934 / HTOCM;
   
    params.num_satellite_temperatures = 11;
    params.satellite_temperatures = (double*) malloc(params.num_satellite_temperatures * sizeof(double));
    params.satellite_temperatures[0]  = 300.0;
    params.satellite_temperatures[1]  = 290.0;
    params.satellite_temperatures[2]  = 280.0;
    params.satellite_temperatures[3]  = 270.0;
    params.satellite_temperatures[4]  = 260.0;
    params.satellite_temperatures[5]  = 250.0;
    params.satellite_temperatures[6]  = 240.0;
    params.satellite_temperatures[7]  = 230.0;
    params.satellite_temperatures[8]  = 220.0;
    params.satellite_temperatures[9]  = 210.0;
    params.satellite_temperatures[10] = 200.0;
    
    params.partial_partition_function_ratios = (double*) malloc(params.num_satellite_temperatures * sizeof(double)); 
    params.partial_partition_function_ratios[0]  = 2.69423e+05; // 300 K
    params.partial_partition_function_ratios[1]  = 2.69477e+05; // 290 K
    params.partial_partition_function_ratios[2]  = 2.69618e+05; // 280 K
    params.partial_partition_function_ratios[3]  = 2.69842e+05; // 270 K
    params.partial_partition_function_ratios[4]  = 2.69944e+05; // 260 K
    params.partial_partition_function_ratios[5]  = 2.70007e+05; // 250 K
    params.partial_partition_function_ratios[6]  = 2.70329e+05; // 240 K
    params.partial_partition_function_ratios[7]  = 2.70396e+05; // 230 K
    params.partial_partition_function_ratios[8]  = 2.70735e+05; // 220 K
    params.partial_partition_function_ratios[9]  = 2.70820e+05; // 210 K
    params.partial_partition_function_ratios[10] = 2.71273e+05; // 200 K
    
    params.cf_filenames = (const char **) malloc(params.num_satellite_temperatures * sizeof(char *));
    params.cf_filenames[0]  = "./CF-CH4-CO2-F-300.0.txt";
    params.cf_filenames[1]  = "./CF-CH4-CO2-F-290.0.txt";
    params.cf_filenames[2]  = "./CF-CH4-CO2-F-280.0.txt";
    params.cf_filenames[3]  = "./CF-CH4-CO2-F-270.0.txt";
    params.cf_filenames[4]  = "./CF-CH4-CO2-F-260.0.txt";
    params.cf_filenames[5]  = "./CF-CH4-CO2-F-250.0.txt";
    params.cf_filenames[6]  = "./CF-CH4-CO2-F-240.0.txt";
    params.cf_filenames[7]  = "./CF-CH4-CO2-F-230.0.txt";
    params.cf_filenames[8]  = "./CF-CH4-CO2-F-220.0.txt";
    params.cf_filenames[9]  = "./CF-CH4-CO2-F-210.0.txt";
    params.cf_filenames[10] = "./CF-CH4-CO2-F-200.0.txt";

    double base_temperature = 300.0;
    CFncArray ca = calculate_correlation_array_and_save(&ms, &params, base_temperature);

    free(params.cf_filenames);
    free(params.satellite_temperatures);
    free(params.partial_partition_function_ratios);

    free_trajectory(&traj);
    free_ms(&ms);
    free_cfnc_array(ca); 

    MPI_Finalize();

    return 0;
}
