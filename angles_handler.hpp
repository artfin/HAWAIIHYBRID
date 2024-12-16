#pragma once

#include <cmath>
#include <vector>

#include <Eigen/Dense>

extern Eigen::Matrix3d SPhi, SPhi_dot;
extern Eigen::Matrix3d STheta, STheta_dot; 
extern Eigen::Matrix3d Sphi1t, Sphi1t_dot;
extern Eigen::Matrix3d Stheta1t, Stheta1t_dot;
extern Eigen::Matrix3d Spsi1t, Spsi1t_dot;
extern Eigen::Matrix3d S1;

extern Eigen::Matrix3d Sphi2t, Sphi2t_dot;
extern Eigen::Matrix3d Stheta2t, Stheta2t_dot;
extern Eigen::Matrix3d Spsi2t, Spsi2t_dot;
extern Eigen::Matrix3d S2;

extern Eigen::Matrix3d Sphiem, Sphiem_dot;
extern Eigen::Matrix3d Sthetaem, Sthetaem_dot;
extern Eigen::Matrix3d Spsiem, Spsiem_dot;

extern Eigen::Vector3d zvec, dd, dd2;
extern Eigen::Vector3d vv1, vv2, vv1dot, vv2dot;

void Sx_filler(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle); 
void Sx_dot_filler(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle);
void Sz_filler(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle); 
void Sz_dot_filler(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle);
void Sy_filler(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle);
void Sy_dot_filler(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle);

void Sy_filler_non_standard(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle);
void Sy_dot_filler_non_standard(Eigen::Ref<Eigen::Matrix3d> m, const double sin_angle, const double cos_angle);

void linear_molecule_atom_lab_to_mol(double *qlab, double *qmol);
void linear_molecule_atom_Jacobi_mol_by_lab(Eigen::Ref<Eigen::MatrixXd> jac, std::vector<double> const& qlab, std::vector<double> const& qmol);

void linear_molecule_linear_molecule_lab_to_mol(std::vector<double> const& qlab, std::vector<double> & qmol);
void linear_molecule_linear_molecule_Jacobi_mol_by_lab(Eigen::Ref<Eigen::MatrixXd> jac, std::vector<double> const& qlab, std::vector<double> const& qmol);

void CH4_linear_molecule_lab_to_kal(std::vector<double> const& qlab, std::vector<double> & qkal);
void CH4_linear_molecule_Jacobi_kal_by_lab(Eigen::Ref<Eigen::MatrixXd> jac, std::vector<double> const& qlab, std::vector<double> const& qkal); 
