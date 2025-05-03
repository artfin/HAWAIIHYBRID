#define USE_MPI
#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

#include "ai_pes_n2_ar_pip_nn.hpp"
#include "ai_ids_n2_ar_pip_nn.hpp"

double pes(double *q)
// PHI THETA R PHI1T THETA1T  
{
    XYZ xyz = lab_to_xyz(q);
    return pes_xyz(&xyz); 
} 

void dpes(double *q, double *dpesdq) 
// PHI THETA R PHI1T THETA1T  
{
    XYZ xyz = lab_to_xyz(q);    
    XYZ dxyz{}; 

    dpes_xyz(&xyz, &dxyz);
    
    double Phi     = q[0];
    double Theta   = q[1];
    double R       = q[2];
    double phi1t   = q[3];
    double theta1t = q[4];
    
    double sinPhi   = sin(Phi);
    double cosPhi   = cos(Phi);
    double sinTheta = sin(Theta);
    double cosTheta = cos(Theta);
    
    /* dVdPhi */   dpesdq[0] = R * (-sinPhi * sinTheta * dxyz.c[6] + cosPhi * sinTheta * dxyz.c[7]) / HTOCM;
    /* dVdTheta */ dpesdq[1] = R * ( cosPhi * cosTheta * dxyz.c[6] + sinPhi * cosTheta * dxyz.c[7] - sinTheta * dxyz.c[8]) / HTOCM;
    /* dVdR */     dpesdq[2] = (cosPhi * sinTheta * dxyz.c[6] + sinPhi * sinTheta * dxyz.c[7] + cosTheta * dxyz.c[8]) / HTOCM;
    
    double sinPhi1t = sin(phi1t);
    double cosPhi1t = cos(phi1t);
    double sinTheta1t = sin(theta1t);
    double cosTheta1t = cos(theta1t);
    
    /* dVdphi1t */   dpesdq[3] = L_N2/2.0 * (-sinPhi1t * sinTheta1t * dxyz.c[0] + cosPhi1t * sinTheta1t * dxyz.c[1] + \
                                              sinPhi1t * sinTheta1t * dxyz.c[3] - cosPhi1t * sinTheta1t * dxyz.c[4]) / HTOCM;
    /* dVdtheta1t */ dpesdq[4] = L_N2/2.0 * (cosPhi1t * cosTheta1t * dxyz.c[0] + sinPhi1t * cosTheta1t * dxyz.c[1] - sinTheta1t * dxyz.c[2] - \
                                             cosPhi1t * cosTheta1t * dxyz.c[3] - sinPhi1t * cosTheta1t * dxyz.c[4] + sinTheta1t * dxyz.c[5]) / HTOCM;
}

void dipole_lab(double *q, double diplab[3]) {
    double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);
    
    Eigen::Vector3d dipmol = dipole_bf(qmol[0], qmol[4]); 
    
    double sinphiem, cosphiem;
    double sinthetaem, costhetaem;
    double sinpsiem, cospsiem;

    sincos(qmol[1], &sinphiem, &cosphiem);
    sincos(qmol[2], &sinthetaem, &costhetaem);
    sincos(qmol[3], &sinpsiem, &cospsiem);

    Sz_filler(Sphiem, sinphiem, cosphiem);
    Sx_filler(Sthetaem, sinthetaem, costhetaem);
    Sz_filler(Spsiem, sinpsiem, cospsiem);
       
    Eigen::Vector3d diplab_eig = Sphiem.transpose() * Sthetaem.transpose() * Spsiem.transpose() * dipmol; 
   
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

    double MU = m_N2 * m_Ar / (m_N2 + m_Ar); 
    double I1[2] = {II_N2, II_N2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL, seeds[_wrank]);
    
    pes_init(true); 
    dipole_init(true);
    dipole = dipole_lab;

    double tolerance = 1e-12;
    Trajectory traj = init_trajectory(&ms, tolerance);
    
    CalcParams params = {};
    params.ps                  = FREE_AND_METASTABLE;
    params.sampler_Rmin        = 4.51;
    params.sampler_Rmax        = 40.0;
    params.niterations         = 100;
    params.total_trajectories  = 10000000;
    params.cvode_tolerance     = 1e-12;
    params.sampling_time       = 200.0;
    params.MaxTrajectoryLength = 65536;
    params.Rcut                = 40.0;
    params.initialM0_npoints   = 30000000;
    params.initialM2_npoints   = 1000;
    params.pesmin              = -97.3 / HTOCM;
   
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
  
    /* 
    params.partial_partition_function_ratios[0]  = 9.60078e+05; // F, 500 K
    params.partial_partition_function_ratios[0]  = 9.13128e+05; // F, 490 K
    params.partial_partition_function_ratios[0]  = 8.67309e+05; // F, 480 K
    params.partial_partition_function_ratios[0]  = 8.22728e+05; // F, 470 K
    params.partial_partition_function_ratios[0]  = 7.79693e+05; // F, 460 K
    params.partial_partition_function_ratios[0]  = 7.38145e+05; // F, 450 K
    params.partial_partition_function_ratios[0]  = 6.97816e+05; // F, 440 K
    params.partial_partition_function_ratios[0]  = 6.58911e+05; // F, 430 K
    params.partial_partition_function_ratios[0]  = 6.21391e+05; // F, 420 K
    params.partial_partition_function_ratios[0]  = 5.84966e+05; // F, 410 K
    params.partial_partition_function_ratios[0]  = 5.50081e+05; // F, 400 K
    params.partial_partition_function_ratios[0]  = 5.16356e+05; // F, 390 K
    params.partial_partition_function_ratios[0]  = 4.83923e+05; // F, 380 K
    params.partial_partition_function_ratios[0]  = 4.52775e+05; // F, 370 K
    params.partial_partition_function_ratios[0]  = 4.22907e+05; // F, 360 K
    params.partial_partition_function_ratios[0]  = 3.94145e+05; // F, 350 K
    params.partial_partition_function_ratios[0]  = 3.66689e+05; // F, 340 K
    params.partial_partition_function_ratios[0]  = 3.40299e+05; // F, 330 K
    params.partial_partition_function_ratios[0]  = 3.15211e+05; // F, 320 K
    params.partial_partition_function_ratios[0]  = 2.91132e+05; // F, 310 K
    */
    params.partial_partition_function_ratios[0]  = 2.68234e+05; // F, 300 K
    params.partial_partition_function_ratios[1]  = 2.46503e+05; // F, 290 K
    params.partial_partition_function_ratios[2]  = 2.25782e+05; // F, 280 K
    params.partial_partition_function_ratios[3]  = 2.06254e+05; // F, 270 K
    params.partial_partition_function_ratios[4]  = 1.87716e+05; // F, 260 K
    params.partial_partition_function_ratios[5]  = 1.70185e+05; // F, 250 K
    params.partial_partition_function_ratios[6]  = 1.53696e+05; // F, 240 K
    params.partial_partition_function_ratios[7]  = 1.38227e+05; // F, 230 K
    params.partial_partition_function_ratios[8]  = 1.23724e+05; // F, 220 K
    params.partial_partition_function_ratios[9]  = 1.10192e+05; // F, 210 K
    params.partial_partition_function_ratios[10] = 9.75779e+04; // F, 200 K
    /*                                                               
    params.partial_partition_function_ratios[11] = 8.58810e+04; // F, 190 K 
    params.partial_partition_function_ratios[11] = 7.50338e+04; // F, 180 K 
    params.partial_partition_function_ratios[12] = 6.50880e+04; // F, 170 K 
    params.partial_partition_function_ratios[13] = 5.59632e+04; // F, 160 K 
    params.partial_partition_function_ratios[14] = 4.76590e+04; // F, 150 K 
    params.partial_partition_function_ratios[15] = 4.01343e+04; // F, 140 K 
    params.partial_partition_function_ratios[16] = 3.33696e+04; // F, 130 K 
    params.partial_partition_function_ratios[17] = 2.73527e+04; // F, 120 K 
    params.partial_partition_function_ratios[18] = 2.20291e+04; // F, 110 K 
    params.partial_partition_function_ratios[19] = 1.73815e+04; // F, 100 K
    params.partial_partition_function_ratios[19] = 1.33802e+04; // F, 90 K
    params.partial_partition_function_ratios[19] = 9.99014e+03; // F, 80 K
    params.partial_partition_function_ratios[19] = 7.17673e+03; // F, 70 K
    */

    /*                                                               
    params.partial_partition_function_ratios[0] = 4.68641e+01; // B, 150 K
    params.partial_partition_function_ratios[0] = 4.76290e+01; // B, 140 K
    params.partial_partition_function_ratios[0] = 4.85049e+01; // B, 130 K
    params.partial_partition_function_ratios[0] = 4.95999e+01; // B, 120 K
    params.partial_partition_function_ratios[0] = 5.09598e+01; // B, 110 K
    params.partial_partition_function_ratios[0] = 5.26148e+01; // B, 100 K
    params.partial_partition_function_ratios[0] = 5.48078e+01; // B, 90 K
    params.partial_partition_function_ratios[0] = 5.77305e+01; // B, 80 K
    params.partial_partition_function_ratios[0] = 6.18422e+01; // B, 70 K
    */

    params.cf_filenames = (const char **) malloc(params.num_satellite_temperatures * sizeof(char *));
    params.cf_filenames[0]  = "./CF-N2-Ar-F-300.0.txt";
    params.cf_filenames[1]  = "./CF-N2-Ar-F-290.0.txt";
    params.cf_filenames[2]  = "./CF-N2-Ar-F-280.0.txt";
    params.cf_filenames[3]  = "./CF-N2-Ar-F-270.0.txt";
    params.cf_filenames[4]  = "./CF-N2-Ar-F-260.0.txt";
    params.cf_filenames[5]  = "./CF-N2-Ar-F-250.0.txt";
    params.cf_filenames[6]  = "./CF-N2-Ar-F-240.0.txt";
    params.cf_filenames[7]  = "./CF-N2-Ar-F-230.0.txt";
    params.cf_filenames[8]  = "./CF-N2-Ar-F-220.0.txt";
    params.cf_filenames[9]  = "./CF-N2-Ar-F-210.0.txt";
    params.cf_filenames[10] = "./CF-N2-Ar-F-200.0.txt";

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
