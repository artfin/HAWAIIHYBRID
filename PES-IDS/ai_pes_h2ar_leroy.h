/*Файл собран согласно статье 
R. J. Le Roy and J. M. Hutson: Potential energy surfaces for H, with Ar, Kr, and Xe, 
J. Chem. Phys., Vol. 86, No.2, 15 January 1987

Input: R - Bohr
Output Enrg - cm-1
Сборка 06.02.2020, QwanXiong
Update 10.10.2021/21.12.2024, artfin
*/
#ifndef AR_H2_PES_LEROY_DERIVATIVES_H_
#define AR_H2_PES_LEROY_DERIVATIVES_H_

#include <math.h>

#include "constants.h"

#ifdef __cplusplus
extern "C" {
#endif

double pes_h2ar(double R, double Theta); 
void dpes_h2ar(double R, double Theta, double *dR, double *dTheta);

double pes_in_cm(double R, double Theta);

double pes_lab(double *q);
void dpes_lab(double *q, double *dpesdq);

#ifdef __cplusplus
}
#endif

#endif
