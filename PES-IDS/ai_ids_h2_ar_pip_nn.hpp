#ifndef AI_IDS_H2_AR_PIP_NN_
#define AI_IDS_H2_AR_PIP_NN_

#include <stdlib.h>
#include <Eigen/Dense>

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

typedef struct {
    double c[15];
} XYZ4q;

void dipole_init(bool log);
Eigen::Vector3d dipole_xyz(XYZ4q *xyz);
Eigen::Vector3d dipole_bf(double R, double Theta); 

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // AI_IDS_H2_AR_PIP_NN_
