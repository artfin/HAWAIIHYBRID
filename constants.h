#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#define Boltzmann  1.380649e-23                                // SI: J * K^(-1)
#define Hartree    4.3597447222071e-18                         // SI: J 
#define HTOCM      2.1947463136320e5                           // 1 Hartree in cm-1
#define HkT        (Hartree/Boltzmann)                         // to use as:  -V[a.u.]*`HkT`/T
#define VkT        (HkT / HTOCM)                               // to use as:  -V[cm-1]*`VkT`/T
    
#define BohrToAng 0.529177210903

#define ZeroCoeff  (0.00361479637/(4.0*M_PI))

#define RAMTOAMU 1822.888485332

// M. Wang, G. Audi, F.G. Kondev, W.J. Huang, S. Naimi, X. Xu. The AME2016 atomic mass evaluation. Tables, graphs and references. 
// http://nuclearmasses.org/resources_folder/Wang_2017_Chinese_Phys_C_41_030003.pdf 
#define m_H  (1.007825032241 * RAMTOAMU)
#define m_C  (12.000000000000 * RAMTOAMU)
#define m_N  (14.003074004460 * RAMTOAMU)
#define m_O  (15.994914619598 * RAMTOAMU)
#define m_Ar (39.9623831237 * RAMTOAMU)

#define m_H2 (2.0 * m_H)
#define m_CO2 (m_C + 2.0 * m_O)

// ? 
#define l_CO2 4.398

// Source: K. P. Huber, G. Herzberg. Molecular Spectra and Molecular Structure. IV. Constants of Diatomic molecules, Springer, US, 1979
// did not found B0 but this value of L_H2 is mentioned in several works, e.g.:
// O. Denis-Alpizar, Y. Kalugina, T. Stoecklin, M. Hernandez Vera, F. lique, A new ab initio potential energy surface for the collisional excitation of HCN by para- and ortho-H2. J. Chem. Phys., 139, 224301, 2013.
// W. Rijks, P. E. S. Wormer, Correlated van der Waals coefficients for dimers consisting of He, Ne, H2, and N2. J. Chem. Phys., 88, 5704, 1988.
#define l_H2 1.448736 
    
#define II_CO2 (m_O/2.0*l_CO2*l_CO2)
#define II_H2  (m_H/2.0*l_H2*l_H2)

#endif // CONSTANTS_H_
