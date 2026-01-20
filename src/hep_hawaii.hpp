#ifndef HAWAII_HEP_HPP_
#define HAWAII_HEP_HPP_

#include "hawaii.h"

// included in <complex.h> and collides with templates in hep
#undef I

#include <hep/mc-mpi.hpp>

typedef double (*Integrand)(hep::mc_point<double> const&);

void transform_variables(hep::mc_point<double> const& x, double* transformed, double* Jac);

double integrand_M0(hep::mc_point<double> const& x);
double integrand_M2(hep::mc_point<double> const& x);
double integrand_pf(hep::mc_point<double> const& x);

void mpi_perform_integration(MoleculeSystem *ms, Integrand integrand, CalcParams *params, double Temperature, size_t niterations, size_t npoints, double *m, double *q);

#endif // HAWAII_HEP_HPP_

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

