#include "constants.h"

#include <math.h>

void dipole_lab(double *q, double diplab[3]) 
// PHI THETA R PHI1T THETA1T  
{
    diplab[0] = mu_CO * sin(q[4]) * cos(q[3]); 
    diplab[1] = mu_CO * sin(q[4]) * sin(q[3]);
    diplab[2] = mu_CO * cos(q[4]);

#if 0
    double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);
    
    double dipmol[3];
    dipmol[0] = arco_dipx_ind(qmol[0] /* R */, qmol[4] /* Theta */);
    dipmol[1] = 0.0; 
    dipmol[2] = arco_dipz_ind(qmol[0] /* R */, qmol[4] /* Theta */);

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
#endif
}
