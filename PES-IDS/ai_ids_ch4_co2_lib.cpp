#include "ai_ids_ch4_co2.hpp"
#include "angles_handler.hpp"

#define NLAB 8
#define NKAL 5

extern "C" {
void dipole_lab(double *qlab, double diplab[3]) {
    double qkal[NKAL];
    CH4_linear_molecule_lab_to_kal(qlab, qkal);
    
    double dipkal[3];
    dipole_vector(qkal, dipkal);

    double sinphi, cosphi;
    double sintheta, costheta;
    double sinpsi, cospsi;

    sincos(qlab[3], &sinphi, &cosphi);
    sincos(qlab[4], &sintheta, &costheta);
    sincos(qlab[5], &sinpsi, &cospsi);

    Sz_filler(Sphi1t, sinphi, cosphi);
    Sx_filler(Stheta1t, sintheta, costheta);
    Sz_filler(Spsi1t, sinpsi, cospsi);
      
    Eigen::Vector3d dipkal_eig = Eigen::Map<Eigen::Vector3d>(dipkal, 3);
    Eigen::Vector3d diplab_eig = Sphi1t.transpose() * Stheta1t.transpose() * Spsi1t.transpose() * dipkal_eig; 
   
    diplab[0] = diplab_eig(0); 
    diplab[1] = diplab_eig(1); 
    diplab[2] = diplab_eig(2); 
}
}
