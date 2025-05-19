#include "c_basis_2_2_1_3_intermolecular.h"

extern "C" void evpoly_2_2_1_3_intermolecular(double* y, Eigen::Ref<Eigen::RowVectorXd> p) {
    p(0) = y[9] + y[8];
    p(1) = y[6] + y[3];
    p(2) = y[8]*y[9];
    p(3) = y[6]*y[9] + y[6]*y[8] + y[3]*y[9] + y[3]*y[8];
    p(4) = y[3]*y[6];
    p(5) = y[9]*y[9] + y[8]*y[8];
    p(6) = y[6]*y[6] + y[3]*y[3];
    p(7) = y[6]*y[8]*y[9] + y[3]*y[8]*y[9];
    p(8) = y[3]*y[6]*y[9] + y[3]*y[6]*y[8];
    p(9) = y[8]*y[9]*y[9] + y[8]*y[8]*y[9];
    p(10) = y[6]*y[9]*y[9] + y[6]*y[8]*y[8] + y[3]*y[9]*y[9] + y[3]*y[8]*y[8];
    p(11) = y[6]*y[6]*y[9] + y[6]*y[6]*y[8] + y[3]*y[3]*y[9] + y[3]*y[3]*y[8];
    p(12) = y[3]*y[6]*y[6] + y[3]*y[3]*y[6];
    p(13) = y[9]*y[9]*y[9] + y[8]*y[8]*y[8];
    p(14) = y[6]*y[6]*y[6] + y[3]*y[3]*y[3];
}
