#include "ai_ids_h2_ar_pip_nn.hpp"

#include "constants.h"

#include "c_basis_2_2_1_3_intermolecular.h"
#include "c_basis_2_1_1_1_3_intermolecular.h"
#include "c_basis_1_1_2_1_3_intermolecular.h"

static const size_t NATOMS = 5;
static const size_t NDIST = NATOMS * (NATOMS - 1) / 2;
static double YIJ[NDIST];

#define QMODEL_IMPLEMENTATION
#include "QModel.hpp"

using namespace qm;

static QModel qmodel;

double sw(double x) {
    double x_i = 6.0;
    double x_f = 20.0;

    if (x < x_i) {
        return 0.0;
    } else if (x < x_f) {
        double r = (x - x_i) / (x_f - x_i);
        double r3 = r * r * r;
        double r4 = r3 * r;
        double r5 = r4 * r;

        return 10.0 * r3 - 15.0 * r4 + 6.0 * r5;
    } else {
        return 1.0;
    }
}

void make_yij_intermolecular_exp5(const double * x, double* yij)
{
    const double a0 = 2.0;
    double drx, dry, drz;
   
    size_t k = 0;
    for (size_t i = 0; i < NATOMS; ++i) {
        for (size_t j = i + 1; j < NATOMS; ++j) {
            
            drx = x[3*i    ] - x[3*j    ]; 
            dry = x[3*i + 1] - x[3*j + 1];
            drz = x[3*i + 2] - x[3*j + 2];
            
            double dst = std::sqrt(drx*drx + dry*dry + drz*drz);
            
            if (i == 0 && j == 1) { yij[k] = 0.0; k = k + 1; continue; } // H1  H2
            if (i == 0 && j == 2) { yij[k] = 0.0; k = k + 1; continue; } // H1  H1'
            if (i == 0 && j == 3) { yij[k] = 0.0; k = k + 1; continue; } // H1  H2'
            if (i == 1 && j == 2) { yij[k] = 0.0; k = k + 1; continue; } // H2  H1'
            if (i == 1 && j == 3) { yij[k] = 0.0; k = k + 1; continue; } // H2  H2'
            if (i == 2 && j == 3) { yij[k] = 0.0; k = k + 1; continue; } // H1' H2'

            double dst5 = dst * dst * dst * dst * dst;
            double s = sw(dst);
            yij[k] = (1.0 - s) * std::exp(-dst / a0) + s * 1e3 / dst5;

            k++;
        }
    }
}

void dipole_init(bool log)
{
    qmodel.init("./PES-IDS/npz/h2-ar-dipoleq-effquad.npz", log);
}

void QMODEL_FILL_POLY(std::vector<size_t> const& npolys, double *x, Eigen::Ref<Eigen::RowVectorXd> p) 
{
    make_yij_intermolecular_exp5(x, YIJ);
    
    evpoly_1_1_2_1_3_intermolecular(YIJ, p.segment(0,                     npolys[0]));
    evpoly_2_1_1_1_3_intermolecular(YIJ, p.segment(npolys[0],             npolys[1]));
    evpoly_2_2_1_3_intermolecular(YIJ,   p.segment(npolys[0] + npolys[1], npolys[2]));
}

void QMODEL_FILL_POLY_INF(std::vector<size_t> const& npolys, Eigen::Ref<Eigen::RowVectorXd> p) {
    p.segment(0, npolys[0])                     = Eigen::RowVectorXd::Zero(npolys[0]);
    p.segment(npolys[0], npolys[1])             = Eigen::RowVectorXd::Zero(npolys[1]);
    p.segment(npolys[0] + npolys[1], npolys[2]) = Eigen::RowVectorXd::Zero(npolys[2]);
}


XYZ4q xyz4q_from_bf(double R, double r, double Theta) 
{
    return XYZ4q {
        .c  = { /* H1  */  r/2.0 * sin(Theta), 0.0,  r/2.0 * cos(Theta),
                /* H2  */ -r/2.0 * sin(Theta), 0.0, -r/2.0 * cos(Theta),
                /* H1' */  r/2.0 * cos(Theta), 0.0, -r/2.0 * sin(Theta),
                /* H2' */ -r/2.0 * cos(Theta), 0.0,  r/2.0 * sin(Theta),
                /* Ar  */  0.0, 0.0, R } 
    };
}

Eigen::Vector3d dipole_bf(double R, double Theta)
{
    XYZ4q xyz = xyz4q_from_bf(R, l_H2, Theta);
    return dipole_xyz(&xyz);
} 

Eigen::Vector3d dipole_xyz(XYZ4q *xyz4q) 
{
    Eigen::RowVectorXd q      = qmodel.forward(xyz4q->c);
    Eigen::RowVectorXd infq   = qmodel.forward_inf();
    Eigen::RowVectorXd q_corr = q - infq;
    Eigen::MatrixXd coords    = Eigen::Map<Eigen::Matrix<double, NATOMS, 3, Eigen::RowMajor>>(xyz4q->c);
        
    return q_corr * coords; 
}

#include "angles_handler.hpp"

void dipole_lab(double *q, double diplab[3]) {
    double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);
    
    Eigen::Vector3d dipmol = dipole_bf(qmol[0], qmol[4]); 
    
    double sinphiem, cosphiem;
    double sinthetaem, costhetaem;
    double sinpsiem, cospsiem;

    sincos(qmol[1], &sinphiem, &cosphiem);
    sincos(qmol[2], &sinthetaem, &costhetaem);
    sincos(qmol[3], &sinpsiem, &cospsiem);

    Sz_filler(Sphiem, sinphiem, cosphiem);
    Sx_filler(Sthetaem, sinthetaem, costhetaem);
    Sz_filler(Spsiem, sinpsiem, cospsiem);
       
    Eigen::Vector3d diplab_eig = Sphiem.transpose() * Sthetaem.transpose() * Spsiem.transpose() * dipmol; 
   
    diplab[0] = diplab_eig(0); 
    diplab[1] = diplab_eig(1); 
    diplab[2] = diplab_eig(2); 
}
