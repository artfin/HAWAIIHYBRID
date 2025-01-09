#include "hawaii.h"

#include "ai_pes_co2ar.h"

#include "ai_ids_co2ar.hpp"
static AI_IDS_co2_ar co2_ar_ids;

#include "angles_handler.hpp"

double pes(double *q) {
    double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);
    return pes_co2ar(qmol[0], qmol[4]);
} 

void dpes(double *q, double *dq) {
    UNUSED(q);
    UNUSED(dq);
    TODO("dpes");
}

#define QP_SIZE 10


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

int main()
{
    uint32_t seed = 43; // mt_goodseed();
    init_pes();
    co2_ar_ids.init();

    dipole = dipole_lab;

    double MU = m_CO2 * m_Ar / (m_CO2 + m_Ar); 
    double I1[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL, seed);

    CalcParams params = {};
    params.ps                               = FREE_AND_METASTABLE;
    params.sampler_Rmin                     = 4.5;
    params.sampler_Rmax                     = 40.0;
    params.initialM0_npoints                = 1000000;
    params.partial_partition_function_ratio = 1.0;
    params.pesmin                           = -195.6337098547 / HTOCM;
    
    double T = 300.0;

    double M0, M0_std;
    calculate_M0(&ms, &params, T, &M0, &M0_std); 

    printf("M0 = %.10e +/- %.10e [%.10e ... %.10e]\n", M0, M0_std, M0-M0_std, M0+M0_std);
    printf("Error: %.3f%%\n", M0_std/M0 * 100.0);

    if (assert_float_is_equal_to(M0, 1.350e-09, 1e-11) > 0) {
        exit(1);
    }


    free_ms(&ms);
    free_pes();

    return 0;
}
