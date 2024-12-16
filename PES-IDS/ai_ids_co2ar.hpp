#pragma once

#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>

#include <gsl/gsl_sf_legendre.h>

struct AI_IDS_co2_ar 
{
public:
	void init();
	void dipole_vector(double R, double Theta, double dip[3]);

    constexpr size_t index(const int l, const int m);
    
    double dipx(double R, double Theta);
    double dipz(double R, double Theta);
    
    const int lmax = 10;
	std::vector<double> LegP;
	
    gsl_sf_legendre_t norm = GSL_SF_LEGENDRE_NONE;
	// to include condon-shortley phase in the definition of GSL associated Legendre polynomial 
	static constexpr double csphase = -1;
    
    bool INITIALIZED = false;
};
