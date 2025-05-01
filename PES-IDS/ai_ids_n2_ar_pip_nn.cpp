#include "ai_ids_n2_ar_pip_nn.hpp"

#include "c_basis_2_2_1_3_purify.h"
#include "c_basis_2_1_1_1_3_purify.h"
#include "c_basis_1_1_2_1_3_purify.h"

#include "constants.h"

static const size_t NATOMS = 5;
static const size_t NDIST = NATOMS * (NATOMS - 1) / 2;
static double YIJ[NDIST];

static double a0 = 2.0;
static double x_i = 9.5;
static double x_f = 12.5;

#define QMODEL_IMPLEMENTATION
#include "QModel.hpp"

using namespace qm;

static QModel qmodel_short;
static QModel qmodel_long;

void make_yij_exp(const double *x, double *yij, double a0)
{
    size_t k = 0;
    for (size_t i = 0; i < NATOMS; ++i) {
        for (size_t j = i + 1; j < NATOMS; ++j) {
            double drx = x[3*i    ] - x[3*j    ]; 
            double dry = x[3*i + 1] - x[3*j + 1];
            double drz = x[3*i + 2] - x[3*j + 2];
            double dst = std::sqrt(drx*drx + dry*dry + drz*drz);

            yij[k] = std::exp(-dst / a0);  
            k++;
        }
    }
}

void QMODEL_FILL_POLY(std::vector<size_t> const& npolys, double *x, Eigen::Ref<Eigen::RowVectorXd> p) 
{
    make_yij_exp(x, YIJ, a0);
    
    evpoly_1_1_2_1_3_purify(YIJ, p.segment(0, npolys[0]));
    evpoly_2_1_1_1_3_purify(YIJ, p.segment(npolys[0], npolys[1]));
    evpoly_2_2_1_3_purify(YIJ, p.segment(npolys[0] + npolys[1], npolys[2]));
}

void QMODEL_FILL_POLY_INF(std::vector<size_t> const& npolys, Eigen::Ref<Eigen::RowVectorXd> p) 
{
    p.segment(0, npolys[0])                     = Eigen::RowVectorXd::Zero(npolys[0]);
    p.segment(npolys[0], npolys[1])             = Eigen::RowVectorXd::Zero(npolys[1]);
    p.segment(npolys[0] + npolys[1], npolys[2]) = Eigen::RowVectorXd::Zero(npolys[2]);
}

double dipole_sw(double x, double x_i, double x_f) 
{
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

void dipole_init(bool log)
{
    qmodel_short.init("./PES-IDS/npz/n2-ar-dipoleq-effquad-short.npz", log);
    qmodel_long.init("./PES-IDS/npz/n2-ar-dipoleq-effquad-long.npz", log);
}

Eigen::Vector3d dipole_from_xyz4q(QModel *qmodel, XYZ4q *xyz4q)
{
    Eigen::RowVectorXd q     = qmodel->forward(xyz4q->c);
    Eigen::RowVectorXd infq  = qmodel->forward_inf();
    Eigen::RowVectorXd qcorr = q - infq;
    Eigen::MatrixXd coords   = Eigen::Map<Eigen::Matrix<double, NATOMS, 3, Eigen::RowMajor>>(xyz4q->c);

    return qcorr * coords;
} 

XYZ4q xyz4q_from_bf(double R, double r, double Theta) 
{
    return XYZ4q {
        .c  = { /* N1  */  r/2.0 * sin(Theta), 0.0,  r/2.0 * cos(Theta),
                /* N2  */ -r/2.0 * sin(Theta), 0.0, -r/2.0 * cos(Theta),
                /* N1' */  r/2.0 * cos(Theta), 0.0, -r/2.0 * sin(Theta),
                /* N2' */ -r/2.0 * cos(Theta), 0.0,  r/2.0 * sin(Theta),
                /* Ar  */  0.0, 0.0, R } 
    };
}

Eigen::Vector3d dipole_bf(double R, double Theta) 
{
    XYZ4q xyz = xyz4q_from_bf(R, L_N2, Theta);
    return dipole_xyz(&xyz);
} 

Eigen::Vector3d dipole_xyz(XYZ4q *xyz)
{
    a0 = 2.0;
    Eigen::Vector3d dipole_short = dipole_from_xyz4q(&qmodel_short, xyz);
    
    a0 = 6.0;
    Eigen::Vector3d dipole_long = dipole_from_xyz4q(&qmodel_long, xyz);

    // we expect here that XYZ4q is prepared using the provided function ???
    // so that it is symmetric and we can safely extract the value of R in the following way

    double nn_com[3];
    nn_com[0] = 0.25 * (xyz->c[0] + xyz->c[3] + xyz->c[6] + xyz->c[9]);
    nn_com[1] = 0.25 * (xyz->c[1] + xyz->c[4] + xyz->c[7] + xyz->c[10]);
    nn_com[2] = 0.25 * (xyz->c[2] + xyz->c[5] + xyz->c[8] + xyz->c[11]);

    double R = (nn_com[0]-xyz->c[12])*(nn_com[0]-xyz->c[12]) + (nn_com[1]-xyz->c[13])*(nn_com[1]-xyz->c[13]) + (nn_com[2]-xyz->c[14])*(nn_com[2]-xyz->c[14]);
    R = sqrt(R);

    double s = dipole_sw(R, x_i, x_f);
    return (1.0 - s) * dipole_short + s * dipole_long;
}

