#include "ai_pes_h2ar_leroy_lib.hpp"
#include "angles_handler.hpp"

#include "ai_pes_h2ar_leroy.h"

extern "C" {
double pes_lab(double *q) {
    static double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);
    return pes_h2ar(qmol[0], qmol[4]);
} 

void dpes_lab(double *q, double *dpesdq) {
    static Eigen::Matrix<double, 5, 5> jac;
    static Eigen::Matrix<double, 5, 1> derivatives_mol, derivatives_lab; 
    static double qmol[5];
    
    jac.setZero();
    linear_molecule_atom_lab_to_mol(q, qmol);
    linear_molecule_atom_Jacobi_mol_by_lab(jac, q, qmol);  
    
    double dR, dTheta;
    dpes_h2ar(qmol[0], qmol[4], &dR, &dTheta); // [R, THETAM] -> [dpes_dR, dpes_dTheta]

    // [dpes_dR, 0, 0, 0, dpes_dTheta]
    derivatives_mol(0) = dR; 
    derivatives_mol(4) = dTheta; 
    
    derivatives_lab = jac * derivatives_mol;
    Eigen::VectorXd::Map(dpesdq, 5) = derivatives_lab;
}
}
