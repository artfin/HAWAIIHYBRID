#define HAWAII_IMPLEMENTATION
#include "hawaii.h"

double pes(double *q) {
    UNUSED(q);
    return 0.0;
}

void dpes(double *q, double *dq) {
    UNUSED(q);
    UNUSED(dq);
    TODO("dpes");
}

#define QP_SIZE 10

int main()
{
    double MU = m_CO2 * m_Ar / (m_CO2 + m_Ar); 
    double I1[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL);

    Array qp = create_array(QP_SIZE);
    double data[] = {7.0, 8.0, 9.0, 10.0, 5.0, 6.0, 11.0, 12.0, 13.0, 14.0};
    init_array(&qp, data, QP_SIZE);

    print_array(qp);

    fill_qp(&ms, qp);

    printf("H = %.10lf\n", Hamiltonian(&ms));

    free_ms(&ms);
    free_array(&qp);

    return 0;
}
