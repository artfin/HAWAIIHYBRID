#pragma once

#include <cmath>
#include <cassert>
#include <vector>

#include "../constants.h"

#include <gsl/gsl_sf_legendre.h>

class AI_PES_co2_ar 
{
public:
    void init();
    double pes(double R, double Theta);
    void dpes(std::vector<double> const& q, std::vector<double> & derivatives);

    // minimum parameters
    const double RMIN = 6.4966051519;
    const double THETAMIN = M_PI / 2.0;
    const double PESMIN = -195.6337098547; // cm-1  

private:
    double dpes_dR    (const double R, const double Theta);
    double dpes_dTheta(const double R, const double Theta);
    
    std::vector<double> legP;
    
    gsl_sf_legendre_t norm = GSL_SF_LEGENDRE_NONE;
	static constexpr double csphase = -1;
    
    bool INITIALIZED = false;
    const int lmax = 10;
};
