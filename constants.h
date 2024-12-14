#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#define Boltzmann  1.380649e-23                                // SI: J * K^(-1)
#define Hartree    4.3597447222071e-18                         // SI: J 
#define HTOCM      2.1947463136320e5                           // 1 Hartree in cm-1
#define HkT        (Hartree/Boltzmann)                         // to use as:  -V[a.u.]*`HkT`/T
#define VkT        (HkT / HTOCM)                               // to use as:  -V[cm-1]*`VkT`/T

#define RAMTOAMU 1822.888485332
    
#define m_C  12.000000000000 * RAMTOAMU
#define m_N  14.003074004460 * RAMTOAMU
#define m_O  15.994914619598 * RAMTOAMU
#define m_Ar 39.9623831237 * RAMTOAMU

#define m_CO2 (m_C + 2.0 * m_O)
#define l_CO2 4.398

#define II_CO2  m_O / 2.0 * l_CO2 * l_CO2

#endif // CONSTANTS_H_
