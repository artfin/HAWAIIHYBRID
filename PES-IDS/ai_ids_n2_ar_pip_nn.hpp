#ifndef AI_IDS_N2_AR_PIP_NN_
#define AI_IDS_N2_AR_PIP_NN_

#include <stdlib.h>
#include <Eigen/Dense>

typedef struct {
    double c[15];
} XYZ4q;

void dipole_init(bool log);
Eigen::Vector3d dipole_xyz(XYZ4q *xyz);
Eigen::Vector3d dipole_bf(double R, double Theta); 

#endif // AI_IDS_N2_AR_PIP_NN_

