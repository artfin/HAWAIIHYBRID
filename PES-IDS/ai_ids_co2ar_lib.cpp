#include "ai_ids_co2ar.hpp"
#include "angles_handler.hpp"

static AI_IDS_co2_ar co2_ar_ids;

extern "C" {
void dipole_init()
{
    co2_ar_ids.init();
}

void dipole_lab(double *q, double diplab[3]) {
    double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);
    
    double dipmol[3];
    co2_ar_ids.dipole_vector(qmol[0], qmol[4], dipmol); 
    
    double sinphiem, cosphiem;
    double sinthetaem, costhetaem;
    double sinpsiem, cospsiem;

    sincos(qmol[1], &sinphiem, &cosphiem);
    sincos(qmol[2], &sinthetaem, &costhetaem);
    sincos(qmol[3], &sinpsiem, &cospsiem);

    Sz_filler(Sphiem, sinphiem, cosphiem);
    Sx_filler(Sthetaem, sinthetaem, costhetaem);
    Sz_filler(Spsiem, sinpsiem, cospsiem);
       
    Eigen::Vector3d dipmol_eig = Eigen::Map<Eigen::Vector3d>(dipmol, 3);
    Eigen::Vector3d diplab_eig = Sphiem.transpose() * Sthetaem.transpose() * Spsiem.transpose() * dipmol_eig; 
   
    diplab[0] = diplab_eig(0); 
    diplab[1] = diplab_eig(1); 
    diplab[2] = diplab_eig(2); 
}
}
