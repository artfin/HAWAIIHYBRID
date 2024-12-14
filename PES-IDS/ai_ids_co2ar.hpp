#pragma once

#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>

#include <Eigen/Dense>
#include <gsl/gsl_sf_legendre.h>

class AI_IDS_co2_ar 
{
public:
	void init();
	Eigen::Vector3d dipole_vector(std::vector<double> const& q);

private:
    constexpr size_t index(const int l, const int m);
    
    double dipx(const double R, const double Theta);
    double dipz(const double R, const double Theta);
    
    const int lmax = 10;
	std::vector<double> LegP;
	
    gsl_sf_legendre_t norm = GSL_SF_LEGENDRE_NONE;
	// to include condon-shortley phase in the definition of GSL associated Legendre polynomial 
	static constexpr double csphase = -1;
    
    bool INITIALIZED = false;
};
