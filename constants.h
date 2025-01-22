#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#define ALU           5.29177210903e-11                           // SI: m
#define Boltzmann     1.380649e-23                                // SI: J * K^(-1)
#define Hartree       4.3597447222071e-18                         // SI: J 
#define HTOCM         2.1947463136320e5                           // 1 Hartree in cm-1
#define ATU           2.4188843265857e-17                         // SI: s 
#define LightSpeed    2.99792458e8                                // SI: m / s
#define LightSpeed_cm 2.99792458e10                               // cm / s
#define ADIPMOMU      8.4783536255e-30                            // SI: C * m
#define EPSILON0      8.8541878128e-12                            // SI: F * m^(-1)
                                            
#define Planck        6.62607015e-34                              // SI: J * s^(-1)
#define HBar          1.054571817e-34                             // SI: J * s^(-1)
#define NL            2.686780111e25                              // SI: m^(-3)

#define HkT        (Hartree/Boltzmann)                         // to use as:  -V[a.u.]*`HkT`/T
#define VkT        (HkT / HTOCM)                               // to use as:  -V[cm-1]*`VkT`/T
    
#define BohrToAng 0.529177210903

#define ZeroCoeff   (0.00361479637/(4.0*M_PI))
#define SecondCoeff 13856114.29344114
    
// MOMENT_SF_COEFF * \int_{-\infty}^{+\infty} J(\nu) * \nu^n d\nu -> cm^(-n-1) amagat^(-2)
#define Moment_SF_Coeff  ((8.0*M_PI*M_PI*M_PI)*NL*NL/3.0/HBar)

#define RAMTOAMU 1822.888485332

// M. Wang, G. Audi, F.G. Kondev, W.J. Huang, S. Naimi, X. Xu. The AME2016 atomic mass evaluation. Tables, graphs and references. 
// http://nuclearmasses.org/resources_folder/Wang_2017_Chinese_Phys_C_41_030003.pdf 
#define m_H  (1.007825032241 * RAMTOAMU)
#define m_C  (12.000000000000 * RAMTOAMU)
#define m_N  (14.003074004460 * RAMTOAMU)
#define m_O  (15.994914619598 * RAMTOAMU)
#define m_Ar (39.9623831237 * RAMTOAMU)

#define m_H2 (2.0 * m_H)
#define m_CO (m_C + m_O)
#define m_CO2 (m_C + 2.0 * m_O)
#define m_CH4 (m_C + 4.0 * m_H)

// ?
#define l_CO 2.132 // bohr 

// ? 
#define l_CO2 4.398

// Source: K. P. Huber, G. Herzberg. Molecular Spectra and Molecular Structure. IV. Constants of Diatomic molecules, Springer, US, 1979
// did not found B0 but this value of L_H2 is mentioned in several works, e.g.:
// O. Denis-Alpizar, Y. Kalugina, T. Stoecklin, M. Hernandez Vera, F. lique, A new ab initio potential energy surface for the collisional excitation of HCN by para- and ortho-H2. J. Chem. Phys., 139, 224301, 2013.
// W. Rijks, P. E. S. Wormer, Correlated van der Waals coefficients for dimers consisting of He, Ne, H2, and N2. J. Chem. Phys., 88, 5704, 1988.
#define l_H2 1.448736 
    
#define II_CO2 (m_O/2.0*l_CO2*l_CO2)
#define II_H2  (m_H/2.0*l_H2*l_H2)
#define II_CO  ((m_C * m_O)/(m_C + m_O)*l_CO*l_CO)

// Source: S. Albert, S. Bauerecker, V. Boudon, L. R. Brown, J. P. Champion, M. Loete, A. Nikitin, M. Quack, Chemical Physics, 2009, 356, 131-146.
// B0 = 5.241040019 cm-1 (page 8; Supplementary material)  
#define l_CH4 2.067354047786849 
// used in our quantum chemistry calculations for CH4-N2 and CH4-CO2
// there is a small inconsistency with B0 from (Albert, 2009)
// l_CH4 -> B0calc = 5.240957524426735
// using this l_CH4 for consistency with PES & IDS 


// source of formula: http://www.pci.tu-bs.de/aggericke/PC4e/exercises/Sol04.pdf 
#define II_CH4 (8.0/3.0 * m_H*l_CH4*l_CH4)
// Product of moments of inertia IIA * IIB * IIC = 33.27803 amu^3 A^6 (amu = unified atomic mass unit = RAM)
// this agrees with the value IIA * IIB * IIC = 33.27341 amu^3 A^6 provided in the
// Computational Chemistry Comparison and Benchmark DataBase (CCCBDB, Release 21, August 2020) by NIST 
// https://cccbdb.nist.gov/exp2x.asp?casno=74828&charge=0 

#endif // CONSTANTS_H_
