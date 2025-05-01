#pragma once

#include <string>
#include <vector>
#include <Eigen/Dense>

#include "scaler.hpp"
#include "cnpy.h"

size_t count_subnets(cnpy::npz_t const& npz);
cnpy::NpyArray get_entry(cnpy::npz_t const& npz, std::string const& entry_name);
std::vector<size_t> parse_architecture_entry(cnpy::npz_t const& npz, std::string const& entry_name);
Eigen::MatrixXd parse_matrixxd_entry(cnpy::npz_t const& npz, std::string const& entry_name); 
Eigen::VectorXd parse_vectorxd_entry(cnpy::npz_t const& npz, std::string const& entry_name);
Eigen::RowVectorXd parse_rowvectorxd_entry(cnpy::npz_t const& npz, std::string const& entry_name);

std::unique_ptr<StandardScaler> build_scaler(cnpy::npz_t& npz, std::string const& scaler_type);
