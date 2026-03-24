#include "angles_handler.hpp"
#include <iostream>

Eigen::Matrix3d SPhi         = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d STheta       = Eigen::Matrix3d::Zero(3, 3);

Eigen::Matrix3d Sphi1t       = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d Stheta1t     = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d Spsi1t       = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d S1           = Eigen::Matrix3d::Zero(3, 3);

Eigen::Matrix3d Sphi2t       = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d Stheta2t     = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d Spsi2t       = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d S2           = Eigen::Matrix3d::Zero(3, 3);

Eigen::Matrix3d SPhi_dot     = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d STheta_dot   = Eigen::Matrix3d::Zero(3, 3);

Eigen::Matrix3d Sphi1t_dot   = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d Stheta1t_dot = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d Spsi1t_dot   = Eigen::Matrix3d::Zero(3, 3);

Eigen::Matrix3d Sphi2t_dot   = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d Stheta2t_dot = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d Spsi2t_dot   = Eigen::Matrix3d::Zero(3, 3);

Eigen::Matrix3d Sphiem       = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d Sthetaem     = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d Spsiem       = Eigen::Matrix3d::Zero(3, 3);

Eigen::Matrix3d Sphiem_dot   = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d Sthetaem_dot = Eigen::Matrix3d::Zero(3, 3);
Eigen::Matrix3d Spsiem_dot   = Eigen::Matrix3d::Zero(3, 3);

Eigen::Vector3d zvec         = Eigen::Vector3d::Unit(3, 2);
Eigen::Vector3d dd           = Eigen::Vector3d::Zero(3, 1);
Eigen::Vector3d dd2          = Eigen::Vector3d::Zero(3, 1);

Eigen::Vector3d vv1          = Eigen::Vector3d::Zero(3, 1);
Eigen::Vector3d vv2          = Eigen::Vector3d::Zero(3, 1);
Eigen::Vector3d vv1dot       = Eigen::Vector3d::Zero(3, 1);
Eigen::Vector3d vv2dot       = Eigen::Vector3d::Zero(3, 1);

void Sx_filler(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle) {
    m(0, 0) = 1.0;
    m(1, 1) = cos_angle; 
    m(1, 2) = sin_angle; 
    m(2, 1) = -sin_angle;
    m(2, 2) = cos_angle; 
}

void Sx_dot_filler(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle) {
    m(1, 1) = -sin_angle;
    m(1, 2) = cos_angle; 
    m(2, 1) = -cos_angle;
    m(2, 2) = -sin_angle;
}

void Sz_filler(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle) {

    m(0, 0) = cos_angle; 
    m(0, 1) = sin_angle; 
    m(1, 0) = -sin_angle;
    m(1, 1) = cos_angle; 
    m(2, 2) = 1.0;
}

void Sz_dot_filler(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle) {
    m(0, 0) = -sin_angle;
    m(0, 1) = cos_angle; 
    m(1, 0) = -cos_angle;
    m(1, 1) = -sin_angle;
}

void Sy_filler(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle) {
    m(0, 0) = cos_angle;
    m(0, 2) = -sin_angle; 
    m(1, 1) = 1.0;
    m(2, 0) = sin_angle;
    m(2, 2) = cos_angle;
}

void Sy_dot_filler(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle) {
    m(0, 0) = -sin_angle;
    m(0, 2) = -cos_angle;
    m(2, 0) = cos_angle;
    m(2, 2) = -sin_angle;
}

void Sy_filler_non_standard(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle) {
    m(0, 0) = cos_angle;
    m(0, 2) = sin_angle; 
    m(1, 1) = 1.0;
    m(2, 0) = -sin_angle;
    m(2, 2) = cos_angle;
}

void Sy_dot_filler_non_standard(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle) {
    m(0, 0) = -sin_angle;
    m(0, 2) = cos_angle;
    m(2, 0) = -cos_angle;
    m(2, 2) = -sin_angle;
}

void linear_molecule_atom_lab_to_mol(double *qlab, double *qmol)
/*
 *  qlab: PHI THETA R PHI1T THETA1T - 
 *  qmol: R PHIEM THETAEM PSIEM THETAM
 */
{
    qmol[0] = qlab[2];
    qmol[1] = qlab[0] + M_PI / 2.0;
    qmol[2] = qlab[1];

    double sinphiem = 0.0, cosphiem = 0.0;
    double sinthetaem = 0.0, costhetaem = 0.0;
    sincos(qmol[1], &sinphiem, &cosphiem);
    sincos(qmol[2], &sinthetaem, &costhetaem);

    double sinphi1t = 0.0, cosphi1t = 0.0;
    double sintheta1t = 0.0, costheta1t = 0.0;
    sincos(qlab[3], &sinphi1t, &cosphi1t);
    sincos(qlab[4], &sintheta1t, &costheta1t);

    Sz_filler(Sphiem, sinphiem, cosphiem);
    Sx_filler(Sthetaem, sinthetaem, costhetaem);

    Sz_filler(Sphi1t, sinphi1t, cosphi1t);
    Sy_filler(Stheta1t, sintheta1t, costheta1t);

    dd = Sthetaem * Sphiem * Sphi1t.transpose() * Stheta1t.transpose() * zvec;

    qmol[4] = std::acos(dd(2));
    double sin_thetam = std::sqrt(1 - dd(2) * dd(2));

    qmol[3] = std::atan2(dd(1) / sin_thetam, dd(0) / sin_thetam);
}

void linear_molecule_atom_Jacobi_mol_by_lab(Eigen::Ref<Eigen::MatrixXd> jac, double *qlab, double *qmol)
/*
 *  qlab: PHI THETA R PHI1T THETA1T -  
 *  qmol: R PHIEM THETAEM PSIEM THETAM
 */
{
    jac(2, 0) = 1.0;
    jac(0, 1) = 1.0;
    jac(1, 2) = 1.0;

    double sinphiem   = 0.0,   cosphiem = 0.0;
    double sinthetaem = 0.0, costhetaem = 0.0;
    double sinpsiem   = 0.0,   cospsiem = 0.0;
    sincos(qmol[1], &sinphiem, &cosphiem);
    sincos(qmol[2], &sinthetaem, &costhetaem);
    sincos(qmol[3], &sinpsiem, &cospsiem);

    double sinphi1t   = 0.0, cosphi1t   = 0.0;
    double sintheta1t = 0.0, costheta1t = 0.0;
    sincos(qlab[3], &sinphi1t, &cosphi1t);
    sincos(qlab[4], &sintheta1t, &costheta1t);

    Sx_dot_filler(Sthetaem_dot, sinthetaem, costhetaem);
    Sz_dot_filler(Sphiem_dot, sinphiem, cosphiem);
    Sz_dot_filler(Sphi1t_dot, sinphi1t, cosphi1t);
    Sy_dot_filler(Stheta1t_dot, sintheta1t, costheta1t);

    double sin_thetam = std::sin(qmol[4]);
  
    //std::cout << "qlab[3]: " << qlab[3] << "\n";
    //std::cout << "qlab[4]: " << qlab[4] << "\n";

    //std::cout << "Sthetaem:\n" << Sthetaem << "\n"; 
    //std::cout << "Sphiem_dot:\n" << Sphiem_dot << "\n"; 
    //std::cout << "Sphi1t:\n" << Sphi1t << "\n"; 
    //std::cout << "Stheta1t:\n" << Stheta1t << "\n"; 
   

    dd = Sthetaem * Sphiem * Sphi1t_dot.transpose() * Stheta1t.transpose() * zvec;
    jac(3, 4) = -dd(2) / sin_thetam; // d(thetam) / d(phi1t)
    jac(3, 3) = (dd(1) * cospsiem - dd(0) * sinpsiem) / sin_thetam; // d(psiem) / d(phi1t)
    
    dd = Sthetaem * Sphiem * Sphi1t.transpose() * Stheta1t_dot.transpose() * zvec;
    jac(4, 4) = -dd(2) / sin_thetam; // d(thetam) / d(theta1t)
    jac(4, 3) = (dd(1) * cospsiem - dd(0) * sinpsiem) / sin_thetam; // d(psiem) / d(theta1t)
    
    dd = Sthetaem * Sphiem_dot * Sphi1t.transpose() * Stheta1t.transpose() * zvec;
    jac(0, 4) = -dd(2) / sin_thetam; // d(thetam) / d(Phi) = d(thetam) / d(phiem)
    jac(0, 3) = (dd(1) * cospsiem - dd(0) * sinpsiem) / sin_thetam; // d(psiem)/d(Phi) = d(psiem)/d(phiem)
    
    dd = Sthetaem_dot * Sphiem * Sphi1t.transpose() * Stheta1t.transpose() * zvec; 
    jac(1, 4) = -dd(2) / sin_thetam; // d(thetam) / d(Theta) = d(thetam) / d(thetaem)
    jac(1, 3) = (dd(1) * cospsiem - dd(0) * sinpsiem) / sin_thetam; // d(psiem)/d(Theta) = d(psiem)/d(thetaem)
}

void linear_molecule_linear_molecule_lab_to_mol(std::vector<double> const& qlab, std::vector<double> & qmol)
/* 
 * qlab:
 *     R PHI THETA
 *     PHI1T THETA1T
 *     PHI2T THETA2T
 * qmol:
 *     R  
 *     PHIEM THETAEM PSIEM 
 *     PHI1M THETA1M (neighbor  molecule)
 *     THETA2M       (reference molecule)
 */
{
    assert(qlab.size() == 7 && "ERROR: expected 7 coordinates");
    assert(qmol.size() == 7 && "ERROR: expected 7 coordinates");
    
    qmol[0] = qlab[0];
    qmol[1] = qlab[1] + M_PI / 2.0;
    qmol[2] = qlab[2];

    double sinphiem, cosphiem;
    double sinthetaem, costhetaem;
    sincos(qmol[1], &sinphiem, &cosphiem);
    sincos(qmol[2], &sinthetaem, &costhetaem);

    double sinphi1t, cosphi1t;
    double sintheta1t, costheta1t;
    sincos(qlab[3], &sinphi1t, &cosphi1t);
    sincos(qlab[4], &sintheta1t, &costheta1t);

    double sinphi2t, cosphi2t;
    double sintheta2t, costheta2t;
    sincos(qlab[5], &sinphi2t, &cosphi2t);
    sincos(qlab[6], &sintheta2t, &costheta2t);

    Sz_filler(Sphiem, sinphiem, cosphiem);
    Sx_filler(Sthetaem, sinthetaem, costhetaem);
    
    Sz_filler(Sphi1t, sinphi1t, cosphi1t);
    Sy_filler(Stheta1t, sintheta1t, costheta1t);

    Sz_filler(Sphi2t, sinphi2t, cosphi2t);
    Sy_filler(Stheta2t, sintheta2t, costheta2t);

    dd = Sthetaem * Sphiem * Sphi2t.transpose() * Stheta2t.transpose() * zvec; 
    
    double sintheta2m = std::sqrt(1 - dd(2) * dd(2));
    qmol[6] = std::acos(dd(2));                                   // theta2m
    qmol[3] = std::atan2(dd(1) / sintheta2m, dd(0) / sintheta2m); // psiem

    double sinpsiem, cospsiem;
    sincos(qmol[3], &sinpsiem, &cospsiem);
    Sz_filler(Spsiem, sinpsiem, cospsiem);

    dd = Spsiem * Sthetaem * Sphiem * Sphi1t.transpose() * Stheta1t.transpose() * zvec;
   
    double sintheta1m = std::sqrt(1.0 - dd(2) * dd(2));
    qmol[5] = std::acos(dd(2));
    qmol[4] = std::atan2(dd(1) / sintheta1m, dd(0) / sintheta1m); 
}

void linear_molecule_linear_molecule_Jacobi_mol_by_lab(Eigen::Ref<Eigen::MatrixXd> jac, std::vector<double> const& qlab, std::vector<double> const& qmol)
/* 
 * qlab:
 *     R PHI THETA
 *     PHI1T THETA1T
 *     PHI2T THETA2T
 * qmol:
 *     R  
 *     PHIEM THETAEM PSIEM 
 *     PHI1M THETA1M (neighbor  molecule)
 *     THETA2M       (reference molecule)
 */
{
    assert(qlab.size() == 7 && "ERROR: expected 7 coordinates");
    assert(qmol.size() == 7 && "ERROR: expected 7 coordinates");
    
    static std::array<Eigen::Matrix3d, 12> mvec;
    static std::array<int, 6> ind;
    
    jac(0, 0) = 1.0; // d(R) / d(R)
    jac(1, 1) = 1.0; // d(phiem) / d(phi)
    jac(2, 2) = 1.0; // d(thetaem) / d(theta)

    double sinphi1m, cosphi1m;
    sincos(qmol[4], &sinphi1m, &cosphi1m); 
    double sintheta1m = std::sin(qmol[5]);
    double sintheta2m = std::sin(qmol[6]);

    double sinphiem, cosphiem;
    double sinthetaem, costhetaem;
    double sinpsiem, cospsiem;
    
    double sinphi1t, cosphi1t;
    double sintheta1t, costheta1t;
    double sinphi2t, cosphi2t;
    double sintheta2t, costheta2t;
    
    sincos(qmol[1], &sinphiem, &cosphiem);
    sincos(qmol[2], &sinthetaem, &costhetaem);
    sincos(qmol[3], &sinpsiem, &cospsiem);
    
    sincos(qlab[3], &sinphi1t, &cosphi1t);
    sincos(qlab[4], &sintheta1t, &costheta1t);
    sincos(qlab[5], &sinphi2t, &cosphi2t);
    sincos(qlab[6], &sintheta2t, &costheta2t);
   
    Sz_filler(Sphiem, sinphiem, cosphiem);
    Sx_filler(Sthetaem, sinthetaem, costhetaem); 
    Sz_filler(Spsiem, sinpsiem, cospsiem);

    Sz_filler(Sphi1t, sinphi1t, cosphi1t);
    Sy_filler(Stheta1t, sintheta1t, costheta1t);
    
    Sz_filler(Sphi2t, sinphi2t, cosphi2t);
    Sy_filler(Stheta2t, sintheta2t, costheta2t);

    Sz_dot_filler(Sphiem_dot, sinphiem, cosphiem);
    Sx_dot_filler(Sthetaem_dot, sinthetaem, costhetaem); 
    Sz_dot_filler(Spsiem_dot, sinpsiem, cospsiem);

    Sz_dot_filler(Sphi1t_dot, sinphi1t, cosphi1t);
    Sy_dot_filler(Stheta1t_dot, sintheta1t, costheta1t);
    
    Sz_dot_filler(Sphi2t_dot, sinphi2t, cosphi2t);
    Sy_dot_filler(Stheta2t_dot, sintheta2t, costheta2t);

    vv1 = Spsiem * Sthetaem * Sphiem * Sphi1t.transpose() * Stheta1t.transpose() * zvec;
    vv2 = Sthetaem * Sphiem * Sphi2t.transpose() * Stheta2t.transpose() * zvec;
    
    mvec = {Sphiem, Sthetaem, Sphi1t, Stheta1t, Sphi2t, Stheta2t, Sphiem_dot, Sthetaem_dot, Sphi1t_dot, Stheta1t_dot, Sphi2t_dot, Stheta2t_dot};
    ind = {0, 0, 0, 0, 0, 0};

    for (size_t k = 0; k < 6; ++k) {
        ind[k] = 6;

        vv2dot = mvec[1 + ind[1]] * mvec[0 + ind[0]] * mvec[4 + ind[4]].transpose() * mvec[5 + ind[5]].transpose() * zvec;
        jac(k + 1, 3) = (vv2dot(1) * cospsiem - vv2dot(0) * sinpsiem) / sintheta2m; // d(psiem) / d(X) 
        jac(k + 1, 6) = -vv2dot(2) / sintheta2m;                                    // d(theta2m) / d(X)

        vv1dot = Spsiem * mvec[1 + ind[1]] * mvec[0 + ind[0]] * mvec[2 + ind[2]].transpose() * mvec[3 + ind[3]].transpose() * zvec; 
        ind[k] = 0;
        vv1dot += jac(k + 1, 3) * Spsiem_dot * mvec[1 + ind[1]] * mvec[0 + ind[0]] * mvec[2 + ind[2]].transpose() * mvec[3 + ind[3]].transpose() * zvec;
	
        jac(k + 1, 4) = (vv1dot(1) * cosphi1m - vv1dot(0) * sinphi1m) / sintheta1m; // d(phi1m) / d(X)
	    jac(k + 1, 5) = -vv1dot(2) / sintheta1m;                                    // d(theta1m) / d(X)
    }

    jac(5, 5) = 0.0;
    jac(6, 5) = 0.0;

    jac(3, 6) = 0.0;
    jac(4, 6) = 0.0;
}


/*
 * This function transforms coordinates in laboratory frame of reference to frame of reference centered on CH4 molecule. It actually does not care how 
 * CH4 is oriented. Initial orientation of CH4 in laboratory frame coincides with the frozen orientation of CH4. 
 * CH4-N2 and CH4-CO2 PES & IDS utilize different orientation of CH4 compared to Y. Kalugina works. So be careful with this. Suffix 'kal' just points to
 * the fact that this frame of reference is centered on CH4. 
 */
void CH4_linear_molecule_lab_to_kal(double *qlab, double *qkal)
/* 
 * qlab:
 *     PHI THETA R
 *     PHI1T THETA1T PSI1T
 *     PHI2T THETA2T
 * qkal:
 *     R PHI1K THETA1K PHI2K THETA2K
 */
{
    qkal[0] = qlab[2];
   
    double sinPhi, cosPhi;
    double sinTheta, cosTheta;
    double sinphi1t, cosphi1t;
    double sintheta1t, costheta1t;
    double sinpsi1t, cospsi1t;
    double sinphi2t, cosphi2t;
    double sintheta2t, costheta2t;

    sincos(qlab[0], &sinPhi, &cosPhi); 
    sincos(qlab[1], &sinTheta, &cosTheta);
    sincos(qlab[3], &sinphi1t, &cosphi1t);
    sincos(qlab[4], &sintheta1t, &costheta1t);
    sincos(qlab[5], &sinpsi1t, &cospsi1t);
    sincos(qlab[6], &sinphi2t, &cosphi2t);
    sincos(qlab[7], &sintheta2t, &costheta2t);
    
    Sz_filler(SPhi, sinPhi, cosPhi);
    Sy_filler_non_standard(STheta, sinTheta, cosTheta);

    Sz_filler(Sphi1t, sinphi1t, cosphi1t);
    Sx_filler(Stheta1t, sintheta1t, costheta1t);
    Sz_filler(Spsi1t, sinpsi1t, cospsi1t);

    Sz_filler(Sphi2t, sinphi2t, cosphi2t);
    Sy_filler_non_standard(Stheta2t, sintheta2t, costheta2t);
    
    S1 = Spsi1t * Stheta1t * Sphi1t;
    dd = S1 * SPhi.transpose() * STheta.transpose() * zvec; 

    // theta1K
    double costheta1k = dd(2);                        
    double sintheta1k = std::sqrt(1 - costheta1k * costheta1k); 
    qkal[2] = std::acos(costheta1k);

    // phi1K 
    double sinphi1k = dd(1) / sintheta1k;
    double cosphi1k = dd(0) / sintheta1k; 
    qkal[1] = std::atan2(sinphi1k, cosphi1k);

    dd = S1 * Sphi2t.transpose() * Stheta2t.transpose() * zvec;

    // theta2K
    double costheta2k = dd(2);                         
    double sintheta2k = std::sqrt(1.0 - dd(2) * dd(2));
    qkal[4] = std::acos(costheta2k);

    // phi2K
    double sinphi2k = dd(1) / sintheta2k; 
    double cosphi2k = dd(0) / sintheta2k; 
    qkal[3] = std::atan2(sinphi2k, cosphi2k); 
}

void CH4_linear_molecule_Jacobi_kal_by_lab(Eigen::Ref<Eigen::MatrixXd> jac, double *qlab, double *qkal)
 /*
  * Матрица упорядочена по столбцам 
  *          d(R)        d(phi1K)    d(theta1K)    d(phi2K)    d(theta2K)
  *  Phi     ...          ...         ...            ...          ...
  *  Theta   ...          ...         ...            ...          ...
  *  R       ...          ...         ...            ...          ...
  *  phi1t   ...          ...         ...            ...          ... 
  *  theta1t ...          ...         ...            ...          ... 
  *  psi1t   ...          ...         ...            ...          ...
  *  phi2t   ...          ...         ...            ...          ...
  *  theta2t ...          ...         ...            ...          ...
  */
{
    jac(2, 0) = 1.0;
    
    
    double sinPhi, cosPhi;
    double sinTheta, cosTheta;
    double sinphi1t, cosphi1t;
    double sintheta1t, costheta1t;
    double sinpsi1t, cospsi1t;
    double sinphi2t, cosphi2t;
    double sintheta2t, costheta2t;

    sincos(qlab[0], &sinPhi, &cosPhi); 
    sincos(qlab[1], &sinTheta, &cosTheta);
    sincos(qlab[3], &sinphi1t, &cosphi1t);
    sincos(qlab[4], &sintheta1t, &costheta1t);
    sincos(qlab[5], &sinpsi1t, &cospsi1t);
    sincos(qlab[6], &sinphi2t, &cosphi2t);
    sincos(qlab[7], &sintheta2t, &costheta2t);
    
    Sz_filler(SPhi, sinPhi, cosPhi);
    Sy_filler_non_standard(STheta, sinTheta, cosTheta);

    Sz_filler(Sphi1t, sinphi1t, cosphi1t);
    Sx_filler(Stheta1t, sintheta1t, costheta1t);
    Sz_filler(Spsi1t, sinpsi1t, cospsi1t);

    Sz_filler(Sphi2t, sinphi2t, cosphi2t);
    Sy_filler_non_standard(Stheta2t, sintheta2t, costheta2t);
    
    Sz_dot_filler(SPhi_dot, sinPhi, cosPhi);
    Sy_dot_filler_non_standard(STheta_dot, sinTheta, cosTheta);

    Sz_dot_filler(Sphi1t_dot, sinphi1t, cosphi1t);
    Sx_dot_filler(Stheta1t_dot, sintheta1t, costheta1t);
    Sz_dot_filler(Spsi1t_dot, sinpsi1t, cospsi1t);

    Sz_dot_filler(Sphi2t_dot, sinphi2t, cosphi2t);
    Sy_dot_filler_non_standard(Stheta2t_dot, sintheta2t, costheta2t);


    std::vector<Eigen::Matrix3d> mvec = {
        SPhi,     STheta,     Sphi1t,     Stheta1t,     Spsi1t,     Sphi2t,     Stheta2t,
        SPhi_dot, STheta_dot, Sphi1t_dot, Stheta1t_dot, Spsi1t_dot, Sphi2t_dot, Stheta2t_dot
    };

    double sinphi1k, cosphi1k;
    double sintheta1k, costheta1k;
    double sinphi2k, cosphi2k;
    double sintheta2k, costheta2k;

    sincos(qkal[1], &sinphi1k, &cosphi1k);
    sincos(qkal[2], &sintheta1k, &costheta1k);
    sincos(qkal[3], &sinphi2k, &cosphi2k);
    sincos(qkal[4], &sintheta2k, &costheta2k);
    
    // jac(IND_LAB, IND_MOL)
    
    // k = 0: d(...)/d(Phi)
    dd  = Spsi1t * Stheta1t * Sphi1t * SPhi_dot.transpose() * STheta.transpose() * zvec;
    dd2 = Spsi1t * Stheta1t * Sphi1t * Sphi2t.transpose() * Stheta2t.transpose() * zvec;
    jac(0, 1) = (-dd(0) * sinphi1k + dd(1) * cosphi1k) / sintheta1k; // d(phi1k)/d(Phi)
    jac(0, 2) = -dd(2) / sintheta1k; // d(theta1k)/d(Phi)
    jac(0, 3) = (-dd2(0) * sinphi2k + dd2(1) * cosphi2k) / sintheta2k; // d(phi2k)/d(Phi)
    jac(0, 4) = -dd2(2) / sintheta2k; // d(theta2k)/d(Phi)
  
    // k = 1: d(...)/d(Theta) 
    dd  = Spsi1t * Stheta1t * Sphi1t * SPhi.transpose() * STheta_dot.transpose() * zvec;
    dd2 = Spsi1t * Stheta1t * Sphi1t * Sphi2t.transpose() * Stheta2t.transpose() * zvec;
    jac(1, 1) = (-dd(0) * sinphi1k + dd(1) * cosphi1k) / sintheta1k; // d(phi1k)/d(Theta)
    jac(1, 2) = -dd(2) / sintheta1k; // d(theta1k)/d(Theta)
    jac(1, 3) = (-dd2(0) * sinphi2k + dd2(1) * cosphi2k) / sintheta2k; // d(phi2k)/d(Theta)
    jac(1, 4) = -dd2(2) / sintheta2k; // d(theta2k)/d(Theta)

    // k = : d(...)/d(phi1t) 
    dd  = Spsi1t * Stheta1t * Sphi1t_dot * SPhi.transpose() * STheta.transpose() * zvec;
    dd2 = Spsi1t * Stheta1t * Sphi1t_dot * Sphi2t.transpose() * Stheta2t.transpose() * zvec;
    jac(3, 1) = (-dd(0) * sinphi1k + dd(1) * cosphi1k) / sintheta1k; // d(phi1k)/d(phi1t)
    jac(3, 2) = -dd(2) / sintheta1k; // d(theta1k)/d(phi1t)
    jac(3, 3) = (-dd2(0) * sinphi2k + dd2(1) * cosphi2k) / sintheta2k; // d(phi2k)/d(phi1t)
    jac(3, 4) = -dd2(2) / sintheta2k; // d(theta2k)/d(phi1t)
    
    // k = 3: d(...)/d(theta1t) 
    dd  = Spsi1t * Stheta1t_dot * Sphi1t * SPhi.transpose() * STheta.transpose() * zvec;
    dd2 = Spsi1t * Stheta1t_dot * Sphi1t * Sphi2t.transpose() * Stheta2t.transpose() * zvec;
    jac(4, 1) = (-dd(0) * sinphi1k + dd(1) * cosphi1k) / sintheta1k; // d(phi1k)/d(theta1t)
    jac(4, 2) = -dd(2) / sintheta1k; // d(theta1k)/d(theta1t)
    jac(4, 3) = (-dd2(0) * sinphi2k + dd2(1) * cosphi2k) / sintheta2k; // d(phi2k)/d(theta1t)
    jac(4, 4) = -dd2(2) / sintheta2k; // d(theta2k)/d(theta1t)
    
    // k = 4: d(...)/d(psi1t) 
    dd  = Spsi1t_dot * Stheta1t * Sphi1t * SPhi.transpose() * STheta.transpose() * zvec;
    dd2 = Spsi1t_dot * Stheta1t * Sphi1t * Sphi2t.transpose() * Stheta2t.transpose() * zvec;
    jac(5, 1) = (-dd(0) * sinphi1k + dd(1) * cosphi1k) / sintheta1k; // d(phi1k)/d(psi1t)
    jac(5, 2) = -dd(2) / sintheta1k; // d(theta1k)/d(psi1t)
    jac(5, 3) = (-dd2(0) * sinphi2k + dd2(1) * cosphi2k) / sintheta2k; // d(phi2k)/d(psi1t)
    jac(5, 4) = -dd2(2) / sintheta2k; // d(theta2k)/d(psi1t)
    
    // k = 5: d(...)/d(phi2t) 
    dd  = Spsi1t * Stheta1t * Sphi1t * SPhi.transpose() * STheta.transpose() * zvec;
    dd2 = Spsi1t * Stheta1t * Sphi1t * Sphi2t_dot.transpose() * Stheta2t.transpose() * zvec;
    jac(6, 1) = (-dd(0) * sinphi1k + dd(1) * cosphi1k) / sintheta1k; // d(phi1k)/d(phi2t)
    jac(6, 2) = -dd(2) / sintheta1k; // d(theta1k)/d(phi2t)
    jac(6, 3) = (-dd2(0) * sinphi2k + dd2(1) * cosphi2k) / sintheta2k; // d(phi2k)/d(phi2t)
    jac(6, 4) = -dd2(2) / sintheta2k; // d(theta2k)/d(phi2t)
    
    // k = 6: d(...)/d(theta2t) 
    dd  = Spsi1t * Stheta1t * Sphi1t * SPhi.transpose() * STheta.transpose() * zvec;
    dd2 = Spsi1t * Stheta1t * Sphi1t * Sphi2t.transpose() * Stheta2t_dot.transpose() * zvec;
    jac(7, 1) = (-dd(0) * sinphi1k + dd(1) * cosphi1k) / sintheta1k; // d(phi1k)/d(theta2t)
    jac(7, 2) = -dd(2) / sintheta1k; // d(theta1k)/d(theta2t)
    jac(7, 3) = (-dd2(0) * sinphi2k + dd2(1) * cosphi2k) / sintheta2k; // d(phi2k)/d(theta2t)
    jac(7, 4) = -dd2(2) / sintheta2k; // d(theta2k)/d(theta2t)

    // theta2k do not depend on Phi, Theta
    jac(0, 4) = 0.0;
    jac(1, 4) = 0.0;

    // theta1k do not depend on phi2t, theta2t
    jac(6, 2) = 0.0;
    jac(7, 2) = 0.0;
} 

void compute_psi_ppsi_for_linear_molecule(double eta, double pEta, double chi, double pChi, double *psi, double *ppsi)
// eta - azimuthal angle
// chi - polar angle
{
    double jx = -pChi * sin(eta) - pEta * cos(eta) / tan(chi);
    double jy = pChi * cos(eta) - pEta * sin(eta) / tan(chi);
    double jz = pEta;

    *ppsi = sqrt(jx*jx + jy*jy + jz*jz); 

    double cosTheta = jz / *ppsi;
    double sinTheta = sqrt(1.0 - jz / *ppsi * jz / *ppsi);
    //double Theta    = acos(cosTheta);

    double Phi    = atan2(jy / *ppsi / sinTheta, jx / *ppsi / sinTheta);
    double sinPhi = sin(Phi);
    double cosPhi = cos(Phi);

    SPhi(0, 0) = -sinPhi;
    SPhi(0, 1) = cosPhi;
    SPhi(0, 2) = 0.0;
    SPhi(1, 0) = -cosPhi;
    SPhi(1, 1) = -sinPhi;
    SPhi(1, 2) = 0.0;
    SPhi(2, 0) = 0.0;
    SPhi(2, 1) = 0.0;
    SPhi(2, 2) = 1.0;
    Sx_filler(STheta, sinTheta, cosTheta);

    dd(0) = cos(eta)*sin(chi);
    dd(1) = sin(eta)*sin(chi);
    dd(2) = cos(chi); 

    dd2 = STheta*SPhi*dd;

    *psi = atan2(dd2(1), dd2(0));
    // printf("psi = %.10e, ppsi = %.10e\n", *psi, *ppsi);

    // save this rotation matrix to rotate vector from molecular to laboratory frame
    S1 = SPhi.transpose() * STheta.transpose(); 
}

void rotate_to_lab_for_linear_molecule(double dipmol[3], double diplab[3])
{
    dd(0) = dipmol[0];
    dd(1) = dipmol[1];
    dd(2) = dipmol[2];

    dd2 = S1 * dd;

    diplab[0] = dd2(0);
    diplab[1] = dd2(1);
    diplab[2] = dd2(2);
}

// H2O-H2O rigid dimer: convert angular coordinates r_ang[9] = {R, Phi, Theta, phi1T, theta1T, psi1T, phi2T, theta2T, psi2T}
// to Cartesian coordinates cart[3][6] (3 coords x 6 atoms: O1 H1 H2 O2 H3 H4)
void h2o_h2o_lab_to_cart(double *r_ang, double cart[3][6])
{
    double R      = r_ang[0];
    double Phi    = r_ang[1];
    double Theta  = r_ang[2];
    double phi1T  = r_ang[3];
    double theta1T = r_ang[4];
    double psi1T  = r_ang[5];
    double phi2T  = r_ang[6];
    double theta2T = r_ang[7];
    double psi2T  = r_ang[8];

    Eigen::Vector3d Rvec(R*sin(Theta)*cos(Phi), R*sin(Theta)*sin(Phi), R*cos(Theta));

    // Equilibrium geometry of a single H2O molecule (in Bohr)
    Eigen::Vector3d Oa( 0.0,              0.0,  0.124885194112977);
    Eigen::Vector3d Ha1( 1.43373587100000, 0.0, -0.991247728887023);
    Eigen::Vector3d Ha2(-1.43373587100000, 0.0, -0.991247728887023);

    Eigen::Vector3d Ob( 0.0,              0.0,  0.124885194112977);
    Eigen::Vector3d Hb1( 1.43373587100000, 0.0, -0.991247728887023);
    Eigen::Vector3d Hb2(-1.43373587100000, 0.0, -0.991247728887023);

    Eigen::Matrix3d Sphi1T_m, Stheta1T_m, Spsi1T_m, Sphi2T_m, Stheta2T_m, Spsi2T_m;

    Sz_filler(Sphi1T_m,  sin(phi1T),  cos(phi1T));
    Sx_filler(Stheta1T_m, sin(theta1T), cos(theta1T));
    Sz_filler(Spsi1T_m,  sin(psi1T),  cos(psi1T));

    Sz_filler(Sphi2T_m,  sin(phi2T),  cos(phi2T));
    Sx_filler(Stheta2T_m, sin(theta2T), cos(theta2T));
    Sz_filler(Spsi2T_m,  sin(psi2T),  cos(psi2T));

    Eigen::Matrix3d S1_m = Sphi1T_m.transpose() * Stheta1T_m.transpose() * Spsi1T_m.transpose();
    Oa  = S1_m * Oa;
    Ha1 = S1_m * Ha1;
    Ha2 = S1_m * Ha2;

    Eigen::Matrix3d S2_m = Sphi2T_m.transpose() * Stheta2T_m.transpose() * Spsi2T_m.transpose();
    Ob  = S2_m * Ob  + Rvec;
    Hb1 = S2_m * Hb1 + Rvec;
    Hb2 = S2_m * Hb2 + Rvec;

    cart[0][0] = Oa(0);  cart[1][0] = Oa(1);  cart[2][0] = Oa(2);
    cart[0][1] = Ha1(0); cart[1][1] = Ha1(1); cart[2][1] = Ha1(2);
    cart[0][2] = Ha2(0); cart[1][2] = Ha2(1); cart[2][2] = Ha2(2);
    cart[0][3] = Ob(0);  cart[1][3] = Ob(1);  cart[2][3] = Ob(2);
    cart[0][4] = Hb1(0); cart[1][4] = Hb1(1); cart[2][4] = Hb1(2);
    cart[0][5] = Hb2(0); cart[1][5] = Hb2(1); cart[2][5] = Hb2(2);
}

// Jacobian d(cart)/d(r_ang) for H2O-H2O rigid dimer
// mat_deriv: 9 x 18 matrix (should be zeroed before call)
void h2o_h2o_der_cart_by_rang(Eigen::Ref<Eigen::MatrixXd> mat_deriv, double cart[3][6], double *r_ang)
{
    double R       = r_ang[0];
    double Phi     = r_ang[1];
    double Theta   = r_ang[2];
    double phi1T   = r_ang[3];
    double theta1T = r_ang[4];
    double psi1T   = r_ang[5];
    double phi2T   = r_ang[6];
    double theta2T = r_ang[7];
    double psi2T   = r_ang[8];

    // Equilibrium geometry of a single H2O molecule (in Bohr)
    Eigen::Vector3d Oa_ini( 0.0,              0.0,  0.124885194112977);
    Eigen::Vector3d Ha1_ini( 1.43373587100000, 0.0, -0.991247728887023);
    Eigen::Vector3d Ha2_ini(-1.43373587100000, 0.0, -0.991247728887023);

    Eigen::Vector3d Ob_ini( 0.0,              0.0,  0.124885194112977);
    Eigen::Vector3d Hb1_ini( 1.43373587100000, 0.0, -0.991247728887023);
    Eigen::Vector3d Hb2_ini(-1.43373587100000, 0.0, -0.991247728887023);

    Eigen::Vector3d Oa, Ha1, Ha2, Ob, Hb1, Hb2, Obd, Hb1d, Hb2d;

    double phi1Tsin = sin(phi1T),   phi1Tcos = cos(phi1T);
    double theta1Tsin = sin(theta1T), theta1Tcos = cos(theta1T);
    double psi1Tsin = sin(psi1T),   psi1Tcos = cos(psi1T);
    double phi2Tsin = sin(phi2T),   phi2Tcos = cos(phi2T);
    double theta2Tsin = sin(theta2T), theta2Tcos = cos(theta2T);
    double psi2Tsin = sin(psi2T),   psi2Tcos = cos(psi2T);

    Eigen::Matrix3d Sphi1T_m, Stheta1T_m, Spsi1T_m, Sphi2T_m, Stheta2T_m, Spsi2T_m;
    Eigen::Matrix3d Sphi1T_dot_m, Stheta1T_dot_m, Spsi1T_dot_m;
    Eigen::Matrix3d Sphi2T_dot_m, Stheta2T_dot_m, Spsi2T_dot_m;

    Sz_filler(Sphi1T_m,       phi1Tsin,   phi1Tcos);
    Sx_filler(Stheta1T_m,     theta1Tsin, theta1Tcos);
    Sz_filler(Spsi1T_m,       psi1Tsin,   psi1Tcos);

    Sz_filler(Sphi2T_m,       phi2Tsin,   phi2Tcos);
    Sx_filler(Stheta2T_m,     theta2Tsin, theta2Tcos);
    Sz_filler(Spsi2T_m,       psi2Tsin,   psi2Tcos);

    Sz_dot_filler(Sphi1T_dot_m,   phi1Tsin,   phi1Tcos);
    Sx_dot_filler(Stheta1T_dot_m, theta1Tsin, theta1Tcos);
    Sz_dot_filler(Spsi1T_dot_m,   psi1Tsin,   psi1Tcos);

    Sz_dot_filler(Sphi2T_dot_m,   phi2Tsin,   phi2Tcos);
    Sx_dot_filler(Stheta2T_dot_m, theta2Tsin, theta2Tcos);
    Sz_dot_filler(Spsi2T_dot_m,   psi2Tsin,   psi2Tcos);

    Eigen::Vector3d Rvec(R*sin(Theta)*cos(Phi), R*sin(Theta)*sin(Phi), R*cos(Theta));
    Eigen::Vector3d Rvec_dR(sin(Theta)*cos(Phi), sin(Theta)*sin(Phi), cos(Theta));
    Eigen::Vector3d Rvec_dT(R*cos(Theta)*cos(Phi), R*cos(Theta)*sin(Phi), -R*sin(Theta));
    Eigen::Vector3d Rvec_dP(-R*sin(Theta)*sin(Phi), R*sin(Theta)*cos(Phi), 0.0);

    // Monomer 2: rotated + translated
    Eigen::Matrix3d S2_m = Sphi2T_m.transpose() * Stheta2T_m.transpose() * Spsi2T_m.transpose();
    Ob  = S2_m * Ob_ini;
    Hb1 = S2_m * Hb1_ini;
    Hb2 = S2_m * Hb2_ini;

    Obd  = Ob  + Rvec;
    Hb1d = Hb1 + Rvec;
    Hb2d = Hb2 + Rvec;

    cart[0][3] = Obd(0);  cart[1][3] = Obd(1);  cart[2][3] = Obd(2);
    cart[0][4] = Hb1d(0); cart[1][4] = Hb1d(1); cart[2][4] = Hb1d(2);
    cart[0][5] = Hb2d(0); cart[1][5] = Hb2d(1); cart[2][5] = Hb2d(2);

    // dR derivative (row 0)
    mat_deriv(0, 9)  = Rvec_dR(0); mat_deriv(0, 10) = Rvec_dR(1); mat_deriv(0, 11) = Rvec_dR(2);
    mat_deriv(0, 12) = Rvec_dR(0); mat_deriv(0, 13) = Rvec_dR(1); mat_deriv(0, 14) = Rvec_dR(2);
    mat_deriv(0, 15) = Rvec_dR(0); mat_deriv(0, 16) = Rvec_dR(1); mat_deriv(0, 17) = Rvec_dR(2);

    // dPhi derivative (row 1)
    mat_deriv(1, 9)  = Rvec_dP(0); mat_deriv(1, 10) = Rvec_dP(1); mat_deriv(1, 11) = Rvec_dP(2);
    mat_deriv(1, 12) = Rvec_dP(0); mat_deriv(1, 13) = Rvec_dP(1); mat_deriv(1, 14) = Rvec_dP(2);
    mat_deriv(1, 15) = Rvec_dP(0); mat_deriv(1, 16) = Rvec_dP(1); mat_deriv(1, 17) = Rvec_dP(2);

    // dTheta derivative (row 2)
    mat_deriv(2, 9)  = Rvec_dT(0); mat_deriv(2, 10) = Rvec_dT(1); mat_deriv(2, 11) = Rvec_dT(2);
    mat_deriv(2, 12) = Rvec_dT(0); mat_deriv(2, 13) = Rvec_dT(1); mat_deriv(2, 14) = Rvec_dT(2);
    mat_deriv(2, 15) = Rvec_dT(0); mat_deriv(2, 16) = Rvec_dT(1); mat_deriv(2, 17) = Rvec_dT(2);

    // Monomer 1: rotated only
    Eigen::Matrix3d S1_m = Sphi1T_m.transpose() * Stheta1T_m.transpose() * Spsi1T_m.transpose();
    Oa  = S1_m * Oa_ini;
    Ha1 = S1_m * Ha1_ini;
    Ha2 = S1_m * Ha2_ini;

    cart[0][0] = Oa(0);  cart[1][0] = Oa(1);  cart[2][0] = Oa(2);
    cart[0][1] = Ha1(0); cart[1][1] = Ha1(1); cart[2][1] = Ha1(2);
    cart[0][2] = Ha2(0); cart[1][2] = Ha2(1); cart[2][2] = Ha2(2);

    // dphi1T derivative (row 3)
    S1_m = Sphi1T_dot_m.transpose() * Stheta1T_m.transpose() * Spsi1T_m.transpose();
    Oa  = S1_m * Oa_ini;
    Ha1 = S1_m * Ha1_ini;
    Ha2 = S1_m * Ha2_ini;

    mat_deriv(3, 0) = Oa(0);  mat_deriv(3, 1) = Oa(1);  mat_deriv(3, 2) = Oa(2);
    mat_deriv(3, 3) = Ha1(0); mat_deriv(3, 4) = Ha1(1); mat_deriv(3, 5) = Ha1(2);
    mat_deriv(3, 6) = Ha2(0); mat_deriv(3, 7) = Ha2(1); mat_deriv(3, 8) = Ha2(2);

    // dtheta1T derivative (row 4)
    S1_m = Sphi1T_m.transpose() * Stheta1T_dot_m.transpose() * Spsi1T_m.transpose();
    Oa  = S1_m * Oa_ini;
    Ha1 = S1_m * Ha1_ini;
    Ha2 = S1_m * Ha2_ini;

    mat_deriv(4, 0) = Oa(0);  mat_deriv(4, 1) = Oa(1);  mat_deriv(4, 2) = Oa(2);
    mat_deriv(4, 3) = Ha1(0); mat_deriv(4, 4) = Ha1(1); mat_deriv(4, 5) = Ha1(2);
    mat_deriv(4, 6) = Ha2(0); mat_deriv(4, 7) = Ha2(1); mat_deriv(4, 8) = Ha2(2);

    // dpsi1T derivative (row 5)
    S1_m = Sphi1T_m.transpose() * Stheta1T_m.transpose() * Spsi1T_dot_m.transpose();
    Oa  = S1_m * Oa_ini;
    Ha1 = S1_m * Ha1_ini;
    Ha2 = S1_m * Ha2_ini;

    mat_deriv(5, 0) = Oa(0);  mat_deriv(5, 1) = Oa(1);  mat_deriv(5, 2) = Oa(2);
    mat_deriv(5, 3) = Ha1(0); mat_deriv(5, 4) = Ha1(1); mat_deriv(5, 5) = Ha1(2);
    mat_deriv(5, 6) = Ha2(0); mat_deriv(5, 7) = Ha2(1); mat_deriv(5, 8) = Ha2(2);

    // dphi2T derivative (row 6)
    S2_m = Sphi2T_dot_m.transpose() * Stheta2T_m.transpose() * Spsi2T_m.transpose();
    Ob  = S2_m * Ob_ini;
    Hb1 = S2_m * Hb1_ini;
    Hb2 = S2_m * Hb2_ini;

    mat_deriv(6, 9)  = Ob(0);  mat_deriv(6, 10) = Ob(1);  mat_deriv(6, 11) = Ob(2);
    mat_deriv(6, 12) = Hb1(0); mat_deriv(6, 13) = Hb1(1); mat_deriv(6, 14) = Hb1(2);
    mat_deriv(6, 15) = Hb2(0); mat_deriv(6, 16) = Hb2(1); mat_deriv(6, 17) = Hb2(2);

    // dtheta2T derivative (row 7)
    S2_m = Sphi2T_m.transpose() * Stheta2T_dot_m.transpose() * Spsi2T_m.transpose();
    Ob  = S2_m * Ob_ini;
    Hb1 = S2_m * Hb1_ini;
    Hb2 = S2_m * Hb2_ini;

    mat_deriv(7, 9)  = Ob(0);  mat_deriv(7, 10) = Ob(1);  mat_deriv(7, 11) = Ob(2);
    mat_deriv(7, 12) = Hb1(0); mat_deriv(7, 13) = Hb1(1); mat_deriv(7, 14) = Hb1(2);
    mat_deriv(7, 15) = Hb2(0); mat_deriv(7, 16) = Hb2(1); mat_deriv(7, 17) = Hb2(2);

    // dpsi2T derivative (row 8)
    S2_m = Sphi2T_m.transpose() * Stheta2T_m.transpose() * Spsi2T_dot_m.transpose();
    Ob  = S2_m * Ob_ini;
    Hb1 = S2_m * Hb1_ini;
    Hb2 = S2_m * Hb2_ini;

    mat_deriv(8, 9)  = Ob(0);  mat_deriv(8, 10) = Ob(1);  mat_deriv(8, 11) = Ob(2);
    mat_deriv(8, 12) = Hb1(0); mat_deriv(8, 13) = Hb1(1); mat_deriv(8, 14) = Hb1(2);
    mat_deriv(8, 15) = Hb2(0); mat_deriv(8, 16) = Hb2(1); mat_deriv(8, 17) = Hb2(2);
}

/*
 *  Copyright (C) 2026 A.Finenko & D.Chistikov
 *  Distributed under the GNU General Public License, version 3
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */       

