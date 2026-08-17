// Shared coordinate reconstruction and NN evaluation for the H2O-Ar wrappers.
// This file evaluates the FULL Liu et al. dipole; the public total and induced
// shared libraries select their semantics in their respective thin wrappers.

#include "PES-IDS/h2o_ar_nn_dms_common.hpp"

#include <utility>

#include <Eigen/Dense>

#include "src/angles_handler.hpp"

extern "C" {
    void dipnn_initx_(void);
    void dipnn_inity_(void);
    void dipnn_initz_(void);

    // xct is Fortran xct(4,3): column-major, so it aliases C double[3][4]
    // indexed [axis][atom], with atoms ordered H1, O, H2, Ar.
    void h2oardipnnx_(double *xct, double *dip, double *dipa, double *dipb, double *dipc);
    void h2oardipnny_(double *xct, double *dip, double *dipa, double *dipb, double *dipc);
    void h2oardipnnz_(double *xct, double *dip, double *dipa, double *dipb, double *dipc);
}

namespace {

// The nets take nine raw Cartesians, so they are neither rotationally covariant
// nor permutationally invariant: they are valid only in the frame used for the fit.
// Reconstruction from the 13462 training geometries gives this convention:
//   - H2O centre of mass at the origin;
//   - Ar on +z;
//   - O in the xz half-plane with x >= 0;
//   - y_A >= 0 for the hydrogen in slot 1;
//   - [(A-O) x (C-O)]_y >= 0 for the hydrogen C in slot 3.
// The last two conditions select the fitted dihedral range [0,90] degrees.

Eigen::Vector3d dms_frame_dipole(const Eigen::Vector3d &H1,
                                 const Eigen::Vector3d &O,
                                 const Eigen::Vector3d &H2,
                                 double R)
{
    double xct[3][4];
    for (int axis = 0; axis < 3; ++axis) {
        xct[axis][0] = H1(axis);
        xct[axis][1] = O(axis);
        xct[axis][2] = H2(axis);
        xct[axis][3] = 0.0;
    }
    xct[2][3] = R;

    // dipa/dipb/dipc are individual ensemble members; only their mean is used.
    double dip[3], dipa, dipb, dipc;
    h2oardipnnx_(&xct[0][0], &dip[0], &dipa, &dipb, &dipc);
    h2oardipnny_(&xct[0][0], &dip[1], &dipa, &dipb, &dipc);
    h2oardipnnz_(&xct[0][0], &dip[2], &dipa, &dipb, &dipc);

    return Eigen::Vector3d(dip[0], dip[1], dip[2]);
}

struct DMSFrame {
    Eigen::Matrix3d U; // maps laboratory vectors to the fitted DMS frame
    Eigen::Vector3d O, H1, H2;
    double R;
    bool reflected;
};

DMSFrame build_dms_frame(double cart[3][4])
{
    const Eigen::Vector3d Olab (cart[0][0], cart[1][0], cart[2][0]);
    const Eigen::Vector3d H1lab(cart[0][1], cart[1][1], cart[2][1]);
    const Eigen::Vector3d H2lab(cart[0][2], cart[1][2], cart[2][2]);
    const Eigen::Vector3d Arlab(cart[0][3], cart[1][3], cart[2][3]);

    DMSFrame frame;
    frame.R = Arlab.norm();

    const Eigen::Vector3d ez = Arlab / frame.R;
    const Eigen::Vector3d ox = Olab - Olab.dot(ez) * ez;

    Eigen::Vector3d ex;
    if (ox.norm() > 1.0e-10) {
        ex = ox.normalized();
    } else {
        // O collinear with Ar: azimuth about the C2 axis is immaterial.
        ex = ez.unitOrthogonal();
    }
    const Eigen::Vector3d ey = ez.cross(ex);

    frame.U.row(0) = ex;
    frame.U.row(1) = ey;
    frame.U.row(2) = ez;

    frame.O  = frame.U * Olab;
    frame.H1 = frame.U * H1lab;
    frame.H2 = frame.U * H2lab;

    // Restore the fitted hydrogen labelling.  Swapping identical H atoms changes
    // the cross-product sign without changing the dipole.  A subsequent y mirror
    // selects y_H1 >= 0; its action on the returned y component is undone below.
    if ((frame.H1 - frame.O).cross(frame.H2 - frame.O)(1) < 0.0) {
        std::swap(frame.H1, frame.H2);
    }

    frame.reflected = false;
    if (frame.H1(1) < 0.0) {
        frame.H1(1) = -frame.H1(1);
        frame.H2(1) = -frame.H2(1);
        frame.O(1)  = -frame.O(1);
        frame.reflected = true;
    }

    return frame;
}

Eigen::Vector3d total_dipole_lab(const DMSFrame &frame)
{
    Eigen::Vector3d dipole =
        dms_frame_dipole(frame.H1, frame.O, frame.H2, frame.R);
    if (frame.reflected) dipole(1) = -dipole(1);
    return frame.U.transpose() * dipole;
}

} // namespace

namespace h2o_ar_nn_dms {

void initialize()
{
    // These routines read weights from PES-IDS/h2o-ar-dms/ relative to the
    // working directory.  HAWAIIHYBRID must therefore run from its repository root.
    dipnn_initx_();
    dipnn_inity_();
    dipnn_initz_();
}

void dipole_lab_at_R(double *q, double R, double diplab[3])
{
    double r_ang[6] = {q[2], q[0], q[1], q[3], q[4], q[5]};
    double cart[3][4]; // [axis][atom], atoms ordered O, H1, H2, Ar
    h2o_ar_lab_to_cart(r_ang, cart);

    DMSFrame frame = build_dms_frame(cart);
    frame.R = R;
    const Eigen::Vector3d dipole = total_dipole_lab(frame);
    diplab[0] = dipole(0);
    diplab[1] = dipole(1);
    diplab[2] = dipole(2);
}

void total_dipole_lab(double *q, double diplab[3])
{
    dipole_lab_at_R(q, q[2], diplab);
}

void embedded_monomer_dipole_lab(double *q, double diplab[3])
{
    dipole_lab_at_R(q, monomer_reference_R_bohr, diplab);
}

} // namespace h2o_ar_nn_dms
