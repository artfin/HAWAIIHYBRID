#include "ai_pes_n2_ar_pip_nn.hpp"

#include "c_basis_2_1_4_purify.h"
#include "c_jac_2_1_4_purify.h"

#include "constants.h"

static int NATOMS = 3;
static int NDIST = NATOMS * (NATOMS - 1) / 2;

static double x_i = 6.0;
static double x_f = 20.0;
static double a0 = 2.0;

double sw(double x, double x_i, double x_f) {
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

double dsw(double x, double x_i, double x_f) {
    if (x < x_i) {
        return 0.0;
    } else if (x < x_f) {
        double t = (x - x_i) / (x_f - x_i);
        double t2 = t * t;
        double t3 = t2 * t;
        double t4 = t3 * t;
        return (30.0 * t2 - 60.0 * t3 + 30.0 * t4) / (x_f - x_i); 
    } else {
        return 0.0;
    }
}

void make_yij_n2_ar_purify(const double * x, double* yij)
{
    double drx, dry, drz;
   
    size_t k = 0;
    for (size_t i = 0; i < (size_t) NATOMS; ++i) {
        for (size_t j = i + 1; j < (size_t) NATOMS; ++j) {
            
            drx = x[3*i    ] - x[3*j    ]; 
            dry = x[3*i + 1] - x[3*j + 1];
            drz = x[3*i + 2] - x[3*j + 2];
            
            double dst = std::sqrt(drx*drx + dry*dry + drz*drz);
            
            if (i == 0 && j == 1) { yij[k] = exp(-dst/2.0); k = k + 1; continue; } // N1 N2 

            double dst6 = dst * dst * dst * dst * dst * dst;
            double s = sw(dst, x_i, x_f);
            yij[k] = (1.0 - s) * std::exp(-dst / 2.0) + s * 1e4 / dst6;

            k++;
        }
    }
}

void make_dydr_n2_ar_purify(Eigen::Ref<Eigen::MatrixXd> dydr, const double* x) 
{
    double drx, dry, drz;
    size_t k = 0;

    for (size_t i = 0; i < (size_t) NATOMS; ++i) {
        for (size_t j = i + 1; j < (size_t) NATOMS; ++j) { 
            drx = x[3*i    ] - x[3*j    ]; 
            dry = x[3*i + 1] - x[3*j + 1];
            drz = x[3*i + 2] - x[3*j + 2];
            
            double dst = std::sqrt(drx*drx + dry*dry + drz*drz);

            if (i == 0 && j == 1) { dydr(k, k) = -1.0/a0 * exp(-dst/2.0); k = k + 1; continue; } // N1 N2
                                                                                   
            double dst6 = dst * dst * dst * dst * dst * dst;
            double dst7 = dst6 * dst;

            double s  = sw(dst, x_i, x_f);
            double ds = dsw(dst, x_i, x_f);
            double ee = std::exp(-dst/a0);

            dydr(k, k) = -ds * ee - (1.0 - s)*ee/a0 + ds * 1e4 / dst6 - s * 6e4 / dst7;
            k++;
        }
    }
}

void EVPOLY(double *y, Eigen::Ref<Eigen::RowVectorXd> p)          { evpoly_2_1_4_purify(y, p); }
void EVPOLY_JAC(Eigen::Ref<Eigen::MatrixXd> jac, double *y)       { evpoly_jac_2_1_4_purify(jac, y); }
void MAKE_YIJ(const double *x, double *y)                         { make_yij_n2_ar_purify(x, y); }
void MAKE_DYDR(Eigen::Ref<Eigen::MatrixXd> dydr, const double *x) { make_dydr_n2_ar_purify(dydr, x); } 

#define MLP_IMPLEMENTATION
#include "mlp.hpp"

static const std::string PES_NPZ = "./PES-IDS/npz/n2-ar-nonrigid-17-32-1-f12.npz"; 
static MLPES model;

XYZ lab_to_xyz(double *q) 
{
    double Phi     = q[0];
    double Theta   = q[1];
    double R       = q[2];
    double phi1t   = q[3];
    double theta1t = q[4];

    XYZ xyz = {
        .c = { 
          /* N1 */  L_N2/2.0*cos(phi1t)*sin(theta1t),  L_N2/2.0*sin(phi1t)*sin(theta1t),  L_N2/2.0*cos(theta1t),
          /* N2 */ -L_N2/2.0*cos(phi1t)*sin(theta1t), -L_N2/2.0*sin(phi1t)*sin(theta1t), -L_N2/2.0*cos(theta1t),
          /* Ar */  R*cos(Phi)*sin(Theta),             R*sin(Phi)*sin(Theta),             R*cos(Theta), 
        }
    };

    return xyz;
}

void pes_init(bool log)
{
    model.init(PES_NPZ, NATOMS, log);
}

double pes_xyz(XYZ *xyz)
{
    assert(PES_NPZ == "./PES-IDS/npz/n2-ar-nonrigid-17-32-1-f12.npz");
    double en_inf = 0.0384490895371528; // this cached value is different for each model
    
    return (model.forward(xyz->c) - en_inf) / HTOCM;
}

void dpes_xyz(XYZ *xyz, XYZ *dxyz)
{
    // NOTE: first we have to call `forward` to cache values
    //       only then we can call `backward` on this configuration 
    model.forward(xyz->c);
    model.backward(xyz->c, dxyz->c);
}
