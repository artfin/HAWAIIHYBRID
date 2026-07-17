#ifndef C_JAC_1_2_1_4_PURIFY_H
#define C_JAC_1_2_1_4_PURIFY_H

#include <Eigen/Dense>
extern "C" void evpoly_jac_1_2_1_4_purify(Eigen::Ref<Eigen::MatrixXf> jac, float* y);

#endif