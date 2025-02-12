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
    
    int seeds[] = {34524213, 20352886, 34249627, 14446993, 65856844, 53218983, 91452731, 43166694, 49914920, 31614602, 44212106, 25700068, 96021096, 41595547, 15145361, 62858278, 75686948, 53370524, 60801159, 93036634, 43643350, 33430677, 51191362, 77677374, 44637994, 96225431, 96728247, 26547561, 42091549, 83904190, 32853682, 56714475, 81442185, 29800713, 83391389, 43932158, 12623570, 17745005, 74524621, 66262234, 80812897, 10189375, 16531916, 36957792, 71300019, 56307328, 16252366, 82968870, 99110762, 23132988, 69426597, 37741864, 82809021, 98537769, 92938416, 85932695, 51181738, 31200502, 78750371, 28980280, 75360893, 27244884, 32760885, 80625563, 57361102, 27955980, 32744304, 30885032, 61644289, 57661732, 32720338, 37280444, 89902840, 28700790, 69382092, 28159346, 47042328, 71604793, 75675316, 56526843, 18983850, 43732077, 25705325, 23544410, 31358496, 15786650, 22945185, 34406691, 34119962, 16242964, 28112522, 97309216, 97766249, 33025954, 11191676, 69854703, 10436926, 33965628, 94552710, 71524683};
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
   
    params.num_satellite_temperatures = 6;
    params.satellite_temperatures = (double*) malloc(params.num_satellite_temperatures * sizeof(double));
    params.satellite_temperatures[0]  = 200.0;
    params.satellite_temperatures[1]  = 190.0;
    params.satellite_temperatures[2]  = 180.0;
    params.satellite_temperatures[3]  = 170.0;
    params.satellite_temperatures[4]  = 160.0;
    params.satellite_temperatures[5]  = 150.0;
    
    params.partial_partition_function_ratios = (double*) malloc(params.num_satellite_temperatures * sizeof(double)); 
    //params.partial_partition_function_ratios[0]  =  2.68736e+05; // 400 K
    //params.partial_partition_function_ratios[1]  =  2.68695e+05; // 390 K
    //params.partial_partition_function_ratios[2]  =  2.68634e+05; // 380 K
    //params.partial_partition_function_ratios[3]  =  2.68936e+05; // 370 K
    //params.partial_partition_function_ratios[4]  =  2.68990e+05; // 360 K
    //params.partial_partition_function_ratios[5]  =  2.69024e+05; // 350 K
    //params.partial_partition_function_ratios[6]  =  2.69049e+05; // 340 K
    //params.partial_partition_function_ratios[7]  =  2.69026e+05; // 330 K
    //params.partial_partition_function_ratios[8]  =  2.69329e+05; // 320 K
    //params.partial_partition_function_ratios[9]  =  2.69332e+05; // 310 K
    //params.partial_partition_function_ratios[0]  = 2.69423e+05; // 300 K
    //params.partial_partition_function_ratios[1]  = 2.69477e+05; // 290 K
    //params.partial_partition_function_ratios[2]  = 2.69618e+05; // 280 K
    //params.partial_partition_function_ratios[3]  = 2.69842e+05; // 270 K
    //params.partial_partition_function_ratios[4]  = 2.69944e+05; // 260 K
    //params.partial_partition_function_ratios[5]  = 2.70007e+05; // 250 K
    //params.partial_partition_function_ratios[6]  = 2.70329e+05; // 240 K
    //params.partial_partition_function_ratios[7]  = 2.70396e+05; // 230 K
    //params.partial_partition_function_ratios[8]  = 2.70735e+05; // 220 K
    //params.partial_partition_function_ratios[9]  = 2.70820e+05; // 210 K
    params.partial_partition_function_ratios[0] = 2.71273e+05; // 200 K
    params.partial_partition_function_ratios[1] = 2.71470e+05; // 190 K
    params.partial_partition_function_ratios[2] = 2.71938e+05; // 180 K
    params.partial_partition_function_ratios[3] = 2.72278e+05; // 170 K
    params.partial_partition_function_ratios[4] = 2.72837e+05; // 160 K
    params.partial_partition_function_ratios[5] = 2.73606e+05; // 150 K
    params.partial_partition_function_ratios[6] = 2.74259e+05; // 140 K
    params.partial_partition_function_ratios[7] = 2.75124e+05; // 130 K
    params.partial_partition_function_ratios[8] = 2.76100e+05; // 120 K
    params.partial_partition_function_ratios[9] = 2.77692e+05; // 110 K
    params.partial_partition_function_ratios[10] = 2.79370e+05; // 100 K
                                                            

    params.cf_filenames = (const char **) malloc(params.num_satellite_temperatures * sizeof(char *));
    params.cf_filenames[0]  = "./CF-CH4-CO2-F-200.0.txt";
    params.cf_filenames[1]  = "./CF-CH4-CO2-F-190.0.txt";
    params.cf_filenames[2]  = "./CF-CH4-CO2-F-180.0.txt";
    params.cf_filenames[3]  = "./CF-CH4-CO2-F-170.0.txt";
    params.cf_filenames[4]  = "./CF-CH4-CO2-F-160.0.txt";
    params.cf_filenames[5]  = "./CF-CH4-CO2-F-150.0.txt";

    double base_temperature = 200.0;
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
