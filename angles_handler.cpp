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

/*
 * This definition of Y-rotation matrix is usually used in quantum mechanical considerations.
 * For example, the basis functions for CH4-N2/CH4-CO2 surfaces utilize this Y-rotation matrix.
 */
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

    double sinphiem = 0.0, cosphiem = 0.0;
    double sinthetaem = 0.0, costhetaem = 0.0;
    double sinpsiem = 0.0, cospsiem = 0.0;
    sincos(qmol[1], &sinphiem, &cosphiem);
    sincos(qmol[2], &sinthetaem, &costhetaem);
    sincos(qmol[3], &sinpsiem, &cospsiem);

    double sinphi1t = 0.0, cosphi1t = 0.0;
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
 * the fact that this frame of reference is just centered on CH4. 
 */
void CH4_linear_molecule_lab_to_kal(std::vector<double> const& qlab, std::vector<double> & qkal)
/* 
 * qlab:
 *     R PHI THETA
 *     PHI1T THETA1T PSI1T
 *     PHI2T THETA2T
 * qkal:
 *     R PHI1K THETA1K PHI2K THETA2K
 */
{
    qkal[0] = qlab[0];
   
    double sinPhi, cosPhi;
    double sinTheta, cosTheta;
    double sinphi1t, cosphi1t;
    double sintheta1t, costheta1t;
    double sinpsi1t, cospsi1t;
    double sinphi2t, cosphi2t;
    double sintheta2t, costheta2t;

    sincos(qlab[1], &sinPhi, &cosPhi); 
    sincos(qlab[2], &sinTheta, &cosTheta);
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

void CH4_linear_molecule_Jacobi_kal_by_lab(Eigen::Ref<Eigen::MatrixXd> jac, std::vector<double> const& qlab, std::vector<double> const& qkal)
 /*
  * Матрица упорядочена по столбцам 
  *          d(R)        d(phi1K)    d(theta1K)    d(phi2K)    d(theta2K)
  *  R       ...          ...         ...            ...          ...
  *  Phi     ...          ...         ...            ...          ...
  *  Theta   ...          ...         ...            ...          ...
  *  phi1t   ...          ...         ...            ...          ... 
  *  theta1t ...          ...         ...            ...          ... 
  *  psi1t   ...          ...         ...            ...          ...
  *  phi2t   ...          ...         ...            ...          ...
  *  theta2t ...          ...         ...            ...          ...
  */
{
    jac(0, 0) = 1.0;
    
    double sinPhi, cosPhi;
    double sinTheta, cosTheta;
    double sinphi1t, cosphi1t;
    double sintheta1t, costheta1t;
    double sinpsi1t, cospsi1t;
    double sinphi2t, cosphi2t;
    double sintheta2t, costheta2t;

    sincos(qlab[1], &sinPhi, &cosPhi); 
    sincos(qlab[2], &sinTheta, &cosTheta);
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

    const int NN = 7;
    int ind[NN] = {0, 0, 0, 0, 0, 0, 0};

    double sinphi1k, cosphi1k;
    double sintheta1k, costheta1k;
    double sinphi2k, cosphi2k;
    double sintheta2k, costheta2k;

    sincos(qkal[1], &sinphi1k, &cosphi1k);
    sincos(qkal[2], &sintheta1k, &costheta1k);
    sincos(qkal[3], &sinphi2k, &cosphi2k);
    sincos(qkal[4], &sintheta2k, &costheta2k);

    for (size_t k = 0; k < NN; ++k) { 
        ind[k] = NN;
        dd = mvec[4 + ind[4]] * mvec[3 + ind[3]] * mvec[2 + ind[2]] * mvec[0 + ind[0]].transpose() * mvec[1 + ind[1]].transpose() * zvec; 
        dd2 = mvec[4 + ind[4]] * mvec[3 + ind[3]] * mvec[2 + ind[2]] * mvec[5 + ind[5]].transpose() * mvec[6 + ind[6]].transpose() * zvec;
        
        jac(k + 1, 1) = (-dd(0) * sinphi1k + dd(1) * cosphi1k) / sintheta1k;
        jac(k + 1, 2) = -dd(2) / sintheta1k;
        jac(k + 1, 3) = (-dd2(0) * sinphi2k + dd2(1) * cosphi2k) / sintheta2k;
        jac(k + 1, 4) = -dd2(2) / sintheta2k;

        ind[k] = 0;
    }

    // why doesn't this happen automatically?

    // phi2k, theta2k do not depend on Phi, Theta
    jac(1, 4) = 0.0;
    jac(2, 4) = 0.0;

    // phi1k, theta1k do not depend on phi2t, theta2t
    jac(6, 2) = 0.0;
    jac(7, 2) = 0.0;
} 



