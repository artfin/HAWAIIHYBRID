#include "parse_cnpy.hpp"

size_t count_subnets(cnpy::npz_t const& npz) {
    size_t counter = 0;
    for (auto it = npz.begin(); it != npz.end(); ++it) {
        std::string name = it->first;
        size_t found = name.find("architecture");
        if (found != std::string::npos) counter++; 
    }

    return counter;
}

cnpy::NpyArray get_entry(cnpy::npz_t const& npz, std::string const& entry_name) {
    for (auto it = npz.begin(); it != npz.end(); ++it) {
        std::string _name = it->first;
        if (_name == entry_name) return it->second; 
    }

    throw std::runtime_error("get entry: could not find entry=" + entry_name); 
}

std::vector<size_t> parse_architecture_entry(cnpy::npz_t const& npz, std::string const& entry_name) {
    cnpy::NpyArray np_subnetwork_arch = get_entry(npz, entry_name);
    assert(np_subnetwork_arch.word_size == sizeof(size_t));

    const size_t* p = np_subnetwork_arch.data<size_t>();
    const size_t sz = np_subnetwork_arch.shape[0];

    return std::vector<size_t>{p, p + sz};
}

Eigen::MatrixXd parse_matrixxd_entry(cnpy::npz_t const& npz, std::string const& entry_name) {
    cnpy::NpyArray np_w = get_entry(npz, entry_name);
    assert(np_w.shape.size() == 2);

    const double* pw = np_w.data<double>();
    Eigen::MatrixXd w = Eigen::MatrixXd(np_w.shape[0], np_w.shape[1]);
    w = Eigen::Map<const Eigen::MatrixXd>(pw, np_w.shape[0], np_w.shape[1]);

    return w;
}

Eigen::RowVectorXd parse_rowvectorxd_entry(cnpy::npz_t const& npz, std::string const& entry_name) {
    cnpy::NpyArray np_b = get_entry(npz, entry_name);
    assert(np_b.shape.size() == 1);

    const double* pb = np_b.data<double>();
    Eigen::RowVectorXd b = Eigen::RowVectorXd(np_b.shape[0]);
    b = Eigen::Map<const Eigen::RowVectorXd>(pb, np_b.shape[0]);

    return b;
}

Eigen::VectorXd parse_vectorxd_entry(cnpy::npz_t const& npz, std::string const& entry_name) {
    cnpy::NpyArray np_vector = get_entry(npz, entry_name);
    assert(np_vector.shape.size() == 1);

    const double* pv = np_vector.data<double>();
    Eigen::VectorXd vector = Eigen::VectorXd(np_vector.shape[0]);
    return Eigen::Map<const Eigen::VectorXd>(pv, np_vector.shape[0]);
}

std::unique_ptr<StandardScaler> build_scaler(cnpy::npz_t& npz, std::string const& scaler_type) {

    assert(scaler_type == "xscaler" || scaler_type == "yscaler");

    auto mean = parse_rowvectorxd_entry(npz, scaler_type + ".mean");
    auto scale = parse_rowvectorxd_entry(npz, scaler_type + ".scale");

    return std::make_unique<StandardScaler>(mean, scale);
}


