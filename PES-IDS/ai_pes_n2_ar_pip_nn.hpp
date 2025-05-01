#ifndef AI_PES_N2_AR_PIP_NN_
#define AI_PES_N2_AR_PIP_NN_

#include <stdlib.h>

typedef struct {
    double c[9];
} XYZ;

XYZ lab_to_xyz(double *q); 

void pes_init(bool log);
double pes_xyz(XYZ *xyz);
void dpes_xyz(XYZ *xyz, XYZ *dxyz);

#endif // AI_PES_N2_AR_PIP_NN_
