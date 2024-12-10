#define HAWAII_IMPLEMENTATION
#include "hawaii.h"

void dpes(double *q, double *dq) {
    UNUSED(q);
    UNUSED(dq);
    TODO("dpes");
}

int main()
{
    double I1[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(LINEAR_MOLECULE, ATOM, I1, NULL);

    free_ms(&ms);

    return 0;
}
