#ifndef HAWAII_HEP_H_
#define HAWAII_HEP_H_

#include "hawaii.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    INTEGRAND_M0,
    INTEGRAND_M2,
    INTEGRAND_PF,
} IntegrandType;

void c_mpi_perform_integration(MoleculeSystem *ms, IntegrandType integrand_type, CalcParams *params, double Temperature, size_t niterations, size_t npoints, double *m, double *q);

#ifdef __cplusplus
}
#endif

#endif // HAWAII_HEP_H_

/*
 *  Copyright (C) 2026 A.Finenko & D.Chistikov 
 *  Distributed under the GNU General Public License, version 3
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */       

