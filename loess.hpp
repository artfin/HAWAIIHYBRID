#ifndef LOESS_H_
#define LOESS_H_

#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "hawaii.h"

/*
 * LOESS: locally weighted polynomial regression
 * Blog post with python implementation:
 * https://towardsdatascience.com/loess-373d43b03564
 * The method is proposed by William S. Cleveland in 1979. 
 * https://www.tandfonline.com/doi/abs/10.1080/01621459.1979.10481038  
 *
 * Apparently, LOESS is equivalent to Savitzky-Golay filtering... 
 *
 * NIST example of LOESS Computation:  
 * https://www.itl.nist.gov/div898/handbook/pmd/section1/dep/dep144.htm 
 */
struct LOESS
{
public:
    typedef std::function<double(double)> weight_func;

    LOESS(double *x, double *y, size_t len);
    
    double estimate(double x, size_t ws, size_t degree);
    std::vector<double> eval_linear_ws(double degree, double xmin, double xmax, size_t npoints, size_t wsmin, double wsstep, size_t wsdelay, double xraw_step);

private:
    std::tuple<double, double> normalize(std::vector<double> & a);
    std::vector<size_t> make_window(std::vector<double> const& distances, size_t window_size);

    void select_by_indices(std::vector<double> const& v, std::vector<size_t> const& indices, std::vector<double> & result);
    void select_by_indices(std::vector<double> const& v, std::vector<size_t> const& indices, Eigen::Ref<Eigen::VectorXd> result);

    std::vector<double> calculate_weights(std::vector<double> const& distances, std::vector<size_t> const& window, weight_func f);
    static double tricubic(double x);

    std::vector<double> m_xx;
    std::vector<double> m_yy;

    double m_minxx, m_maxxx;
    double m_minyy, m_maxyy;
};
    
#endif // LOESS_H_
