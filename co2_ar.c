#define HAWAII_IMPLEMENTATION
#include "hawaii.h"

#include "ai_pes_co2ar.hpp"
static AI_PES_co2_ar co2_ar_pes;

// #include "ai_ids_co2ar.hpp"
// static AI_PES_co2_ar ids;

#include "angles_handler.hpp"

double pes(double *q) {
    double qmol[5];
    AnglesHandling::linear_molecule_atom_lab_to_mol(q, qmol);
    return co2_ar_pes.pes(qmol[0], qmol[4]);
} 

void dpes(double *q, double *dq) {
    UNUSED(q);
    UNUSED(dq);
    TODO("dpes");
}

#define QP_SIZE 10

double pesmin = 1.0;

void dipole(double *q, double dip[3]) {
    dip[0] = q[0];
    dip[1] = q[1];
    dip[2] = q[2];
}

int main()
{
    uint32_t seed = mt_goodseed();
    co2_ar_pes.init();
    // ids.init();

    double MU = m_CO2 * m_Ar / (m_CO2 + m_Ar); 
    double I1[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL, seed);

    Array qp = create_array(QP_SIZE);
    double data[] = {7.0, 8.0, 9.0, 10.0, 5.0, 6.0, 11.0, 12.0, 13.0, 14.0};
    init_array(&qp, data, QP_SIZE);

    print_array(qp);

    fill_qp(&ms, qp);

    printf("H = %.10lf\n", Hamiltonian(&ms));

    CalcParams params = {};
    params.ps = BOUND;
    params.sampler_Rmin = 4.5;
    params.sampler_Rmax = 40.0;
    params.initialM0_npoints = 10000;
    
    double T = 300.0;
    double result = calculate_M0(&ms, &params, T); 

    printf("result = %.10lf\n", result);


    free_ms(&ms);
    free_array(&qp);

    // q_generator(MoleculeSystem *ms, CalcParams *params) 

    return 0;
}
