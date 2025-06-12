#include "angles_handler.hpp"
#include "constants.h"

// В работе Pedersen написано, что потенциал предполагает, что theta_m = 0 отвечает C-O-Ar, а theta_m = Pi  —  O-C-Ar, где theta_m определен как угол между вектором от C к О и вектором, 
// соединяющим центр масс CO и Ar. В лабораторной системе координат мы вправе двумя способами выбрать углы (phi1t, theta1t). В ситуации theta1t = 0 молекула ориентирована вдоль оси OZ и 
// положительную координату имеет или  атом C или атом O: (+)C — O(-) или (-) С — O (+). Пусть мы считаем, что theta1t = 0 отвечает ситуации (-) C — O (+). 
// То есть phi1t = (произвольное значение), theta1t = 0, Phi = (произвольное значение), Theta = 0 обозначает конфигурацию Ar-O-C, где все атомы лежат на оси OZ лабораторной системы координат. 
// Применяя пересчет углов в подвижную систему, получаем theta_m = 0, что согласуется с конвенцией, принятой в потенциале Pederson'a.

// Dipole Rizzo
// θ = 0 deg corresponds to Ar-O-C an θ = 180 deg corresponds Ar-C-O
// Dipole direction: -CO+

extern "C" {
	// IPOT = 3 : Use the parameters of the CCSD(T)/aug-cc-pVQZ-33211
	//            3D surface at R1. This surface is probably not
	//            trustworthy beyond the interval (this is NOT tested!)
	//               1.898 bohr <= R1 <= 2.234 bohr
	//            due to the single-reference nature of CCSD(T).
	//            Thus, intended only for 2-D calculations at
	//            a given CO distance in the above trust interval.

    // R1    - CO distance in bohr
    // R2    - distance from CO center-of-mass to Ar in bohr
    // XCOS2 - cos to the intermolecular angle
	void potv(double* res, double* r1, double* r2, double* xcos2);
	void potv_d(double* v, double* vd, double* r1, double* r1d, double* r2, double* r2d,
		        double* xcos2, double* xcos2d);
}

extern "C" {
double pes_lab(double *q) {
    static double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);

    // r1 - CO bond length in bohr
    // r2 - distance from CO center-of-mass to Ar in bohr
    // theta - angle in rad
	double _r1 = l_CO;
	double _r2 = qmol[0];
	double _xcos2 = cos(qmol[4]);

    double r;
	potv(&r, &_r1, &_r2, &_xcos2);
	return r;
} 

void dpes_lab(double *q, double *dpesdq) {
    static Eigen::Matrix<double, 5, 5> jac;
    static Eigen::Matrix<double, 5, 1> derivatives_mol, derivatives_lab; 
    static double qmol[5];
     
    jac.setZero();
    linear_molecule_atom_lab_to_mol(q, qmol);
    linear_molecule_atom_Jacobi_mol_by_lab(jac, q, qmol);  
    
	double _r1 = l_CO; 
	double _r1d = 0.0;
	double _r2 = qmol[0]; 
	double _r2d = 1.0;
	double _xcos2 = cos(qmol[4]);
	double _xcos2d = 0.0;

	double t, dR, dTheta;
	
    potv_d(&t, &dR, &_r1, &_r1d, &_r2, &_r2d, &_xcos2, &_xcos2d);

	_r1d = 0.0;
	_r2d = 0.0;
	_xcos2d = 1.0;
	potv_d(&t, &dTheta, &_r1, &_r1d, &_r2, &_r2d, &_xcos2, &_xcos2d);

    // [dpes_dR, 0, 0, 0, dpes_dTheta]
    derivatives_mol(0) = dR; 
    derivatives_mol(4) = -sin(qmol[4])*dTheta; 

    //for (size_t i = 0; i < 5; ++i) {
    //    printf("derivatives_mol(%zu) = %.10lf\n", i, derivatives_mol(i));
    //}
    //
    //for (size_t i = 0; i < 5; ++i) {
    //    for (size_t j = 0; j < 5; ++j) {
    //        printf("%.5e ", jac(i, j)); 
    //    }
    //    printf("\n");
    //}
    
    derivatives_lab = jac * derivatives_mol;
    Eigen::VectorXd::Map(dpesdq, 5) = derivatives_lab;
}
}

