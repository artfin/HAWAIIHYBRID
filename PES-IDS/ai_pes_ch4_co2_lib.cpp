#include "ai_pes_ch4_co2.h"
#include "angles_handler.hpp"

#define NLAB 8
#define NKAL 5

extern "C" {
double pes_lab(double *qlab) {
    static double qkal[NKAL];
    CH4_linear_molecule_lab_to_kal(qlab, qkal); 
    return pes_ch4co2(qkal[0], qkal[1], qkal[2], qkal[3], qkal[4]); 
}

void dpes_lab(double *qlab, double *dpesdq) {
    static Eigen::Matrix<double, NLAB, NKAL> jac;
    static Eigen::Matrix<double, NKAL, 1> derivatives_kal;
    static Eigen::Matrix<double, NLAB, 1> derivatives_lab; 
    static double qkal[NKAL];

    jac.setZero();
    CH4_linear_molecule_lab_to_kal(qlab, qkal);
    CH4_linear_molecule_Jacobi_kal_by_lab(jac, qlab, qkal);  

    double dR, dphi1, dphi2, dtheta1, dtheta2; 
    dpes_ch4co2(qkal[0], qkal[1], qkal[2], qkal[3], qkal[4], &dR, &dphi1, &dtheta1, &dphi2, &dtheta2); 

    derivatives_kal(0) = dR; 
    derivatives_kal(1) = dphi1;
    derivatives_kal(2) = dtheta1;
    derivatives_kal(3) = dphi2;
    derivatives_kal(4) = dtheta2;

    derivatives_lab = jac * derivatives_kal;

    Eigen::VectorXd::Map(dpesdq, NLAB) = derivatives_lab;
}
}
