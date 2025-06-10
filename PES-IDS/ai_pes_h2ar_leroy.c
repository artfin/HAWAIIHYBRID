#include "ai_pes_h2ar_leroy.h"

#include "angles_handler.hpp"

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


double pes_h2ar(double R, double Theta) {
    return pes_in_cm(R, Theta) / HTOCM;
} 

double dpes_dR(double Rbohr, double Theta);
double dpes_dTheta(double Rbohr, double Theta);

void dpes_h2ar(double R, double Theta, double *dR, double *dTheta)
{
    *dR     = dpes_dR(R, Theta) / HTOCM;
    *dTheta = dpes_dTheta(R, Theta) / HTOCM;
}

double pes_in_cm(double Rbohr, double Theta)
{
	double t1, t4, t9, t12, t15, t18, t21, t24, t27, t28, t31, t33, t52, t56;

    double R = Rbohr * BohrToAng;
	t1 = R * R;
	t4 = exp(-0.35815000000000000e1 * R);
	t9 = t1 * R;
	t12 = t1 * t1;
	t15 = t12 * R;
	t18 = t12 * t1;
	t21 = t12 * t9;
	t24 = t12 * t12;
	t27 = cos(Theta);
	t28 = t27 * t27;
	t31 = t28 * t4;
	t33 = -0.2555000000000000e16 * t1 + 0.74665070545000000e17 * t4 * R + 0.13626147507845875e18 * t4 * t1 + 0.16877397933116667e18 * t4 * t9 + 0.15930933885583086e18 * t4 * t12 + 0.12193846258704485e18 * t4 * t15 + 0.78625828439263148e17 * t4 * t18 + 0.43813155564488664e17 * t4 * t21 + 0.35492851688729762e18 * t4 * t24 - 0.405000000000000e15 * t28 * t1 + 0.9598530000000000e16 * t31;
	t52 = -0.20847430000000000e17 + 0.13400094514989060e18 * t31 * t24 + 0.61965854850446250e17 * t31 * t1 + 0.74943908048957748e17 * t31 * t9 + 0.68401649822148042e17 * t31 * t12 + 0.50236492836750768e17 * t31 * t15 + 0.30912512477084406e17 * t31 * t18 + 0.16384404515581361e17 * t31 * t21 + 0.34377135195000000e17 * t31 * R + 0.20847430000000000e17 * t4 - 0.9598530000000000e16 * t28;
	t56 = 0.50000000000000000e-10 * (t33 + t52) / t24;

    return t56;
}

double dpes_dR(double Rbohr, double Theta)
{
	double t1, t3, t4, t9, t14, t17, t20, t23, t26, t29, t32, t34, t52, t58, t63;

    double R = Rbohr * BohrToAng;
	t1 = R * R;
	t3 = cos(Theta);
	t4 = t3 * t3;
	t9 = exp(-0.35815000000000000e1 * R);
	t14 = t1 * R;
	t17 = t1 * t1;
	t20 = t17 * R;
	t23 = t17 * t1;
	t26 = t17 * t14;
	t29 = t17 * t17;
	t32 = t4 * t9;
	t34 = -0.7665000000000000e16 * t1 - 0.1215000000000000e16 * t4 * t1 - 0.38394120000000000e17 * t4 + 0.29866028218000000e18 * t9 * R + 0.54249090031383500e18 * t9 * t1 + 0.66594518482466668e18 * t9 * t14 + 0.62085068119894844e18 * t9 * t17 + 0.46819089243664639e18 * t9 * t20 + 0.29698713031701371e18 * t9 * t23 + 0.16270578005985482e18 * t9 * t26 + 0.78458408327108100e17 * t9 * t29 + 0.38394120000000000e17 * t32;
	t52 = t29 * R;
	t58 = -0.83389720000000000e17 + 0.83389720000000000e17 * t9 + 0.29340372386277300e17 * t32 * t29 + 0.24745841940178500e18 * t32 * t1 + 0.29832512469583099e18 * t32 * t14 + 0.27100910298296717e18 * t32 * t17 + 0.19784499367413776e18 * t32 * t20 + 0.12087351202449584e18 * t32 * t23 + 0.63548783976129579e17 * t32 * t26 + 0.13750854078000000e18 * t32 * R + 0.23996219252716659e18 * t52 * t4 * t9 + 0.63558824161592820e18 * t52 * t9;
	t63 = -0.10000000000000000e-9 * (t34 + t58) / t52;

    return t63 * BohrToAng;
}

double dpes_dTheta(double Rbohr, double Theta)
{
	double t1, t2, t3, t5, t7, t10, t16, t32, t36;

    double R = Rbohr * BohrToAng;
	t1 = R * R;
	t2 = t1 * t1;
	t3 = t2 * t2;
	t5 = cos(Theta);
	t7 = sin(Theta);
	t10 = exp(-0.35815000000000000e1 * R);
	t16 = t1 * R;
	t32 = -0.405000000000000e15 * t1 + 0.9598530000000000e16 * t10 + 0.13400094514989060e18 * t10 * t3 + 0.61965854850446250e17 * t10 * t1 + 0.74943908048957750e17 * t10 * t16 + 0.68401649822148040e17 * t10 * t2 + 0.50236492836750770e17 * t10 * t2 * R + 0.30912512477084406e17 * t10 * t2 * t1 + 0.16384404515581361e17 * t10 * t2 * t16 + 0.34377135195000000e17 * t10 * R - 0.9598530000000000e16;
	t36 = -0.10000000000000000e-9 / t3 * t5 * t7 * t32;

    return t36;
}

