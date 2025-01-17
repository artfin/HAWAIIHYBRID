#include "loess.hpp"

LOESS::LOESS(double *x, double *y, size_t len) {
    m_xx = std::vector<double>(x, x + len);
    m_yy = std::vector<double>(y, y + len);

    std::tie(m_minxx, m_maxxx) = normalize(m_xx);
    std::tie(m_minyy, m_maxyy) = normalize(m_yy); 
}

std::tuple<double, double> LOESS::normalize(std::vector<double> & a) {

    size_t size = a.size();
    double min_value = *std::min_element(a.begin(), a.end());
    double max_value = *std::max_element(a.begin(), a.end());

    for (size_t k = 0; k < size; ++k) {
        a[k] = (a[k] - min_value) / (max_value - min_value);
    } 
    
    return std::make_tuple(min_value, max_value);
}

std::vector<size_t> LOESS::make_window(std::vector<double> const& distances, size_t window_size) {

    size_t n_distances = distances.size();
    size_t min_idx = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));

    if (min_idx == 0) {
        double *ls = linspace(0, window_size - 1, window_size);
        return std::vector<size_t>(ls, ls + window_size);
    } else if (min_idx == n_distances - 1) {
        double *ls = linspace(n_distances - window_size, n_distances - 1, window_size);
        return std::vector<size_t>(ls, ls + window_size);
    }

    std::vector<size_t> window;
    window.push_back(min_idx);
    
    size_t ind_upper, ind_lower;

    while (window.size() < window_size) {
        // current lower and upper bounds of the window 
        ind_lower = window[0];
        ind_upper = window.back();

        if (ind_lower == 0) {
            window.push_back(ind_upper + 1);
        } else if (ind_upper == n_distances - 1) {
            window.insert(window.begin(), ind_lower - 1);
        } else if (distances[ind_lower - 1] < distances[ind_upper + 1]) {
            window.insert(window.begin(), ind_lower - 1);
        } else {
            window.push_back(ind_upper + 1);
        }
    }
    
    return window;
}

std::vector<double> LOESS::calculate_weights(std::vector<double> const& distances, std::vector<size_t> const& window, weight_func f) {

    size_t window_size = window.size();
    std::vector<double> weights(window_size);

    std::vector<double> window_distances(window_size); 
    select_by_indices(distances, window, window_distances);

    double max_distance = *std::max_element(window_distances.begin(), window_distances.end());

    for (size_t k = 0; k < window_size; ++k) {
        weights[k] = f(window_distances[k] / max_distance);
    }

    return weights;
}

double LOESS::tricubic(double x) {
    return (1.0 - x * x * x) * (1.0 - x * x * x) * (1.0 - x * x * x);  
} 

void LOESS::select_by_indices(std::vector<double> const& v, std::vector<size_t> const& indices, std::vector<double> & result) {
    
    assert(result.size() == indices.size() && "ERROR: expected vectors of equal size");

    size_t size = indices.size();
    for (size_t k = 0; k < size; ++k) {
        result[k] = v[indices[k]];
    }
}

void LOESS::select_by_indices(std::vector<double> const& v, std::vector<size_t> const& indices, Eigen::Ref<Eigen::VectorXd> result) {
    
    assert(static_cast<size_t>(result.size()) == indices.size() && "ERROR: expected vectors of equal size");

    size_t size = indices.size();
    for (size_t k = 0; k < size; ++k) {
        result(k) = v[indices[k]];
    }
}

std::vector<double> LOESS::eval_linear_ws(double degree, double xmin, double xmax, size_t npoints, size_t wsmin, double wsstep, size_t wsdelay, double xraw_step)
{
    std::vector<double> sm(npoints);

    double xstep = (xmax - xmin) / (npoints - 1);
    size_t ws = wsmin;

    for (size_t k = 0; k < npoints; ++k) {
        double x = xmin + xstep * k;
       
        if (k > wsdelay) ws = static_cast<size_t>(wsmin + wsstep*k); 
        assert(ws >= 1); 

        sm[k] = estimate(x, ws, degree);

        if (k % 100 == 0) {
            std::cout << "  k=" << k << "; freq=" << x << "; val = " << sm[k] << "; [npoints] window_size=" << ws << "; [freq] window_size=" << ws * xraw_step << "\n";
        }
    } 

    return sm;
}

double LOESS::estimate(double x, size_t ws, size_t degree)
/*
 * ws [window_size]: number of points to be included in the local averaging 
 */ 
{
    double n_x = (x - m_minxx) / (m_maxxx - m_minxx);

    std::vector<double> distances(m_xx.size());
    for (size_t i = 0; i < distances.size(); ++i) {
        distances[i] = std::abs(m_xx[i] - n_x);
    } 
   
    std::vector<size_t> window = make_window(distances, ws); 
    size_t wss = window.size();
    
    auto weights = calculate_weights(distances, window, tricubic);

    Eigen::VectorXd weights_ = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());
    Eigen::MatrixXd wm = weights_.asDiagonal();
    Eigen::MatrixXd xm = Eigen::MatrixXd::Constant(wss, degree + 1, 1.0);
    
    Eigen::VectorXd m_xx_ = Eigen::VectorXd::Zero(wss, 1);
    Eigen::VectorXd m_yy_ = Eigen::VectorXd::Zero(wss, 1);

    select_by_indices(m_xx, window, m_xx_);
    select_by_indices(m_yy, window, m_yy_);
    
    for (size_t k = 1; k < degree + 1; ++k) {
        xm.col(k) = m_xx_.array().pow(k);
    }

    auto xmt_wm = xm.transpose() * wm;
    auto beta = (xmt_wm * xm).completeOrthogonalDecomposition().pseudoInverse() * xmt_wm * m_yy_;

    Eigen::VectorXd xp = Eigen::VectorXd::Zero(degree + 1);
    for (size_t p = 0; p < degree + 1; ++p) {
        xp(p) = std::pow(n_x, p);
    }

    double y = beta.dot(xp);
    
    return y * (m_maxyy - m_minyy) + m_minyy; 
}

