#include "hep_hawaii.hpp"

static MoleculeSystem* gms = NULL;
static CalcParams *gparams = NULL;
static double gT = 0.0;

void transform_variables(hep::mc_point<double> const& x, double* qp, double* jac)
/*
 * Input: n-dimensional point in [0; 1]^n  hypercube
 * Output: [std::vector<double>] point in phase space
 *         [double]              cumulative jacobian of the transform
 */
{
    assert(x.point().size() == gms->QP_SIZE);

    *jac = 1;

    qp[IR] = 1.0 / x.point()[0]; 
    *jac *= qp[IR] * qp[IR];  
    
    qp[IPR] =  std::tan(M_PI * (x.point()[IPR] - 0.5));
    *jac *= M_PI * (1.0 + qp[IPR] * qp[IPR]);
    
    qp[IPHI] = x.point()[IPHI] * 2.0 * M_PI;
    *jac *= 2.0 * M_PI;
    
    qp[IPPHI] =  std::tan(M_PI * (x.point()[IPPHI] - 0.5));
    *jac *= M_PI * (1.0 + qp[IPPHI] * qp[IPPHI]);
    
    qp[ITHETA] = x.point()[ITHETA] * M_PI;
    *jac *= M_PI;
    
    qp[IPTHETA] =  std::tan(M_PI * (x.point()[IPTHETA] - 0.5));
    *jac *= M_PI * (1.0 + qp[IPTHETA] * qp[IPTHETA]);

    switch (gms->m1.t) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
          qp[6 + IPHI] = x.point()[6 + IPHI] * 2.0 * M_PI;
          *jac *= 2.0 * M_PI;
        
          qp[6 + IPPHI] = std::tan(M_PI * (x.point()[6 + IPPHI] - 0.5));
          *jac *= M_PI * (1.0 + qp[6 + IPPHI] * qp[6 + IPPHI]);
          
          qp[6 + ITHETA] = x.point()[6 + ITHETA] * 2.0 * M_PI;
          *jac *= M_PI;
        
          qp[6 + IPTHETA] = std::tan(M_PI * (x.point()[6 + IPTHETA] - 0.5));
          *jac *= M_PI * (1.0 + qp[6 + IPTHETA] * qp[6 + IPTHETA]);
          break;
        }
        default: {
          TODO("transform_variables");
        }
    }

    switch (gms->m2.t) {
        case ATOM: break; 
        default: {
          TODO("transform_variables");
        }
    }
}

double integrand_pf(hep::mc_point<double> const& x)
{
    assert(gms != NULL);
    assert(gparams != NULL);
    assert(gT > 0);
    
    double jac;
    double qp[gms->QP_SIZE];
    transform_variables(x, qp, &jac);
    put_qp_into_ms(gms, (Array) {.data = qp, .n = gms->QP_SIZE});
    
    double R = qp[IR];
    if (R < gparams->sampler_Rmin || R > gparams->sampler_Rmax) {
        return 0.0;
    }
    
    double energy = Hamiltonian(gms); 
    
    if (gparams->ps == FREE_AND_METASTABLE) {
        if (energy < 0.0) return 0.0; 
    }

    if (gparams->ps == BOUND) {
        if (energy > 0.0) return 0.0; 
    }
    
    return jac * std::exp(-energy * HkT / gT);
}

double integrand_M0(hep::mc_point<double> const& x) 
{
    assert(gms != NULL);
    assert(gparams != NULL);
    assert(gT > 0);

    double jac;
    double qp[gms->QP_SIZE];
    transform_variables(x, qp, &jac);
    put_qp_into_ms(gms, (Array) {.data = qp, .n = gms->QP_SIZE});

    double R = qp[IR];
    if (R < gparams->sampler_Rmin || R > gparams->sampler_Rmax) {
        return 0.0;
    }

    double energy = Hamiltonian(gms); 
    
    if (gparams->ps == FREE_AND_METASTABLE) {
        if (energy < 0.0) return 0.0; 
    }

    if (gparams->ps == BOUND) {
        if (energy > 0.0) return 0.0; 
    }
    
    double dip[3];
    extract_q_and_write_into_ms(gms);
    (*dipole)(gms->intermediate_q, dip);
    
    double dipsq = dip[0]*dip[0] + dip[1]*dip[1] + dip[2]*dip[2]; 

    //std::cout << "R: " << R << " => jac = " << jac << ", energy = " << energy << "\n";

    return jac * dipsq * std::exp(-energy * HkT / gT);
}

void mpi_perform_integration(MoleculeSystem *ms, Integrand integrand, CalcParams *params, double T, size_t niterations, size_t npoints, double *m, double *q)
{
    gms = ms;
    gparams = params;
    gT = T;

    auto results = hep::mpi_vegas(
        MPI_COMM_WORLD,
        hep::make_integrand<double>(integrand, ms->QP_SIZE),
        std::vector<size_t>(niterations, npoints)
    );
    
    hep::mc_result<double> result = hep::accumulate<hep::weighted_with_variance>(results.begin() + 1, results.end());    
    // double chi_square_dof = hep::chi_square_dof<hep::weighted_with_variance>(results.begin() + 1, results.end());

    *m = result.value();
    *q = result.error();
}


