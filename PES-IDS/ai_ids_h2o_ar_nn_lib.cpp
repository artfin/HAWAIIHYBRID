// Shared library wrapper for the H2O-Ar NN dipole moment surface
// Exposes dipole_lab, dipole_init for HAWAIIHYBRID
//
// Model: the three-component NN of PES-IDS/h2o-ar-dms (dipx.f, dipy.f, dipz.f), from
//   Liu, Wang, Zhou, Xie, "A Full-Dimensional ab initio Intermolecular Potential Energy
//   Surface and Dipole Moment Surfaces for H2O-Ar", Current Chinese Science 2, 325 (2022),
//   doi:10.2174/2210298102666220404103308.
// The DMSs are MP2/AVTZ, fitted by NN over an R grid of 4.0-20.0 a0 with a flexible H2O
// (Radau r1, r2 in [1.4, 2.5] a0, Radau angle th1 in [7, 173] deg). Each component is an
// ensemble of three 9 -> 50 -> 50 -> 1 tansig nets whose outputs are averaged. Values are
// in atomic units of dipole, which is what HAWAIIHYBRID expects, so no conversion happens
// here. The paper quotes fitting errors of 6.192 mD (X) and 6.509 mD (Z); it reports
// nothing for Y, and the authors shipped reference data for X and Z only.
//
// INDUCED vs TOTAL: the surface fits the TOTAL dipole of the complex, which at
// large R tends to the H2O permanent dipole (0.727 a.u. = 1.85 D). CIA needs the
// induced part only, so we return
//     mu_ind = mu_tot(R) - mu_tot(R_MONOMER_REF)
// evaluated at one and the same H2O orientation. The subtrahend is the free
// monomer dipole as this fit represents it: it depends on the H2O orientation but
// not on R, and taking it from the net itself (rather than an analytic H2O dipole)
// cancels the fit's own systematic error in the monomer limit.
//
// RANGE: the fit covers R in [4.0, 20] bohr and is not damped or cut off outside it.
// mu_ind vanishes by construction at R = R_MONOMER_REF; beyond it the nets extrapolate
// and the difference grows spuriously (~3e-3 a.u. already at R = 22). This is intended
// for bound states, which rarely sample R > 20.
//
// ACCURACY: mu_ind is a difference of two net evaluations, so the fit error enters it
// twice and does not fully cancel. Measured against the paper's own numbers: at the
// global minimum this wrapper returns 1.8529 D where the paper reports 1.853 D. But at
// R = R_MONOMER_REF the nets give 1.858-1.899 D for the free monomer (true value 1.853 D
// at every orientation), so the subtrahend carries roughly 5-45 mD of error that lands
// directly in mu_ind as an orientation-dependent systematic. Against ~300 mD of induced
// dipole at contact that is a few per cent, growing as mu_ind decays.
//
// Coordinate mapping:
//   HAWAIIHYBRID intermediate_q = [Phi, Theta, R, phi1T, theta1T, psi1T]
//   PES code r_ang              = [R,   Phi,   Theta, phi1T, theta1T, psi1T]
// Monomer 1 is the H2O rotor, monomer 2 is the Ar atom.

#include <cstdio>
#include <utility>

#include <Eigen/Dense>

#include "src/angles_handler.hpp"

// The nets take nine raw Cartesians, so they are neither rotationally covariant
// nor permutationally invariant: they are only valid in the frame the fit was made
// in, and they return the dipole in that same frame. Reconstructed from all 13462
// training geometries, that frame is:
//   - H2O centre of mass at the origin (max offset over the training set: 2.4e-7)
//   - Ar on the +z axis
//   - O in the xz half-plane with x >= 0
//   - y_A >= 0, where A is the hydrogen in slot 1 (holds for all 13462 points)
//   - m_y >= 0, where m = (A - O) x (C - O) and C is the hydrogen in slot 3
//     (holds for all 13462 points; min over the set is -2.2e-8, i.e. roundoff)
//
// The last two together are what confine the fit to the dihedral range phi in
// [0, 90] deg that it was trained on: writing y_A ~ sin(phi) and m_y ~ cos(phi),
// y_A >= 0 gives phi in [0,180] and m_y >= 0 gives phi in [-90,90]. Enforcing only
// the y_A rule leaves phi unconstrained in [0,180]; the nets then extrapolate badly
// and return up to 4 D for a monomer whose dipole is 1.85 D.
static const double R_MONOMER_REF = 20.0; // bohr; upper edge of the fitted R range

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

// Evaluate the total dipole, in the DMS frame, for an H2O already expressed in that
// frame and an Ar placed at (0, 0, R).
static Eigen::Vector3d dms_frame_dipole(const Eigen::Vector3d &H1, const Eigen::Vector3d &O,
                                        const Eigen::Vector3d &H2, double R)
{
    double xct[3][4];
    for (int a = 0; a < 3; ++a) {
        xct[a][0] = H1(a);
        xct[a][1] = O(a);
        xct[a][2] = H2(a);
        xct[a][3] = 0.0;
    }
    xct[2][3] = R;

    // dipa/dipb/dipc are the individual ensemble members; only their mean is used.
    double dip[3], dipa, dipb, dipc;
    h2oardipnnx_(&xct[0][0], &dip[0], &dipa, &dipb, &dipc);
    h2oardipnny_(&xct[0][0], &dip[1], &dipa, &dipb, &dipc);
    h2oardipnnz_(&xct[0][0], &dip[2], &dipa, &dipb, &dipc);

    return Eigen::Vector3d(dip[0], dip[1], dip[2]);
}

// The lab geometry expressed in the fitted frame: U maps lab -> DMS, Od/Ha/Hb are the
// H2O atoms in that frame with the fitted labelling applied, and `reflected` records
// whether the mirror image was taken (in which case dy must be negated on the way out).
struct DMSFrame {
    Eigen::Matrix3d U;
    Eigen::Vector3d Od, Ha, Hb;
    double R;
    bool reflected;
};

static DMSFrame build_dms_frame(double cart[3][4])
{
    Eigen::Vector3d Olab (cart[0][0], cart[1][0], cart[2][0]);
    Eigen::Vector3d H1lab(cart[0][1], cart[1][1], cart[2][1]);
    Eigen::Vector3d H2lab(cart[0][2], cart[1][2], cart[2][2]);
    Eigen::Vector3d Arlab(cart[0][3], cart[1][3], cart[2][3]);

    DMSFrame f;

    // h2o_ar_lab_to_cart already places the H2O centre of mass at the origin and
    // measures Ar from it, so the frames differ by a rotation only.
    f.R = Arlab.norm();

    Eigen::Vector3d ez = Arlab / f.R;
    Eigen::Vector3d ox = Olab - Olab.dot(ez) * ez;

    Eigen::Vector3d ex;
    double oxn = ox.norm();
    if (oxn > 1e-10) {
        ex = ox / oxn;
    } else {
        // O collinear with Ar: the C2 axis lies along R, the azimuth is undefined and
        // the dipole does not depend on it, so any ex perpendicular to ez will do.
        ex = ez.unitOrthogonal();
    }
    Eigen::Vector3d ey = ez.cross(ex);

    // Rows of U are the DMS axes in lab coordinates, so U maps lab -> DMS. The triad
    // is right-handed by construction (ex x ey = ez), hence U is a proper rotation.
    f.U.row(0) = ex;
    f.U.row(1) = ey;
    f.U.row(2) = ez;

    f.Od = f.U * Olab;
    f.Ha = f.U * H1lab;
    f.Hb = f.U * H2lab;

    // Restore the fitted labelling (y_A >= 0 and m_y >= 0). Two operations suffice,
    // and each has a known effect on the pair of signs:
    //   swap A and C   : flips m_y and flips y_A. The hydrogens are identical
    //                    particles, so this is pure relabelling and the dipole is
    //                    unchanged.
    //   reflect y -> -y: flips y_A and leaves m_y alone (for the improper map
    //                    M = diag(1,-1,1), (Ma) x (Mb) = -M(a x b)). This sends the
    //                    geometry to its mirror image, whose dipole is the mirrored
    //                    dipole, so dy has to be undone afterwards.
    // Applying the swap first and the reflection second lands in phi in [0,90], since
    // the reflection cannot disturb the m_y sign the swap just fixed.
    if ((f.Ha - f.Od).cross(f.Hb - f.Od)(1) < 0.0) std::swap(f.Ha, f.Hb);

    f.reflected = false;
    if (f.Ha(1) < 0.0) {
        f.Ha(1) = -f.Ha(1);
        f.Hb(1) = -f.Hb(1);
        f.Od(1) = -f.Od(1); // zero by construction, mirrored for completeness
        f.reflected = true;
    }

    return f;
}

// Total dipole in the lab frame, with Ar placed at R_eval along its actual direction.
// R_eval = f.R gives mu_tot(R); R_eval = R_MONOMER_REF gives the free-monomer dipole
// at the same H2O orientation.
static Eigen::Vector3d dipole_total_lab(const DMSFrame &f, double R_eval)
{
    Eigen::Vector3d d = dms_frame_dipole(f.Ha, f.Od, f.Hb, R_eval);
    if (f.reflected) d(1) = -d(1);
    return f.U.transpose() * d;
}

extern "C" {

void dipole_init(bool log)
{
    // Each of these reads its weights from PES-IDS/h2o-ar-dms/ using a path relative
    // to the working directory, so the driver must run from the repository root.
    dipnn_initx_();
    dipnn_inity_();
    dipnn_initz_();

    if (log) {
        printf("H2O-Ar NN dipole surface initialized (induced dipole, "
               "monomer reference at R = %.1f bohr)\n", R_MONOMER_REF);
    }
}

void dipole_lab(double *q, double diplab[3])
{
    double r_ang[6] = {q[2], q[0], q[1], q[3], q[4], q[5]};

    double cart[3][4]; // [axis][atom], atoms ordered O, H1, H2, Ar
    h2o_ar_lab_to_cart(r_ang, cart);

    DMSFrame f = build_dms_frame(cart);

    Eigen::Vector3d dip_ind = dipole_total_lab(f, f.R) - dipole_total_lab(f, R_MONOMER_REF);

    diplab[0] = dip_ind(0);
    diplab[1] = dip_ind(1);
    diplab[2] = dip_ind(2);
}

} // extern "C"
