#ifndef MY_MATRIX_CALC
#define MY_MATRIX_CALC

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <chrono>
#include <map>
#include <utility>

void read_baypass_file(const std::string file_name, std::vector<int>& total_count, 
  std::vector<std::vector<double>>& data_cov_mat, double maf_filter, std::ofstream& log);
void estimate_pop_cov_mat(const std::vector<int>& total_count, 
  const std::vector<std::vector<double>>& sample_cov_mat, 
  std::vector<std::vector<double>>& pop_cov_mat, std::ofstream& log);


void calculate_equilibrium_moment(const std::vector<std::vector<double>>& mig_mat, 
  Eigen::VectorXd& x, Eigen::SparseMatrix<double>& A_mat, const std::vector<int>& index_order);
void calculate_covariance_structure(const int pop_num, const Eigen::VectorXd& x, 
  std::vector<std::vector<double>>& pop_cov);
void calculate_covariance_sample_structure(const std::vector<std::vector<double>>& pop_cov, 
  const std::vector<int>& sampled_num, std::vector<std::vector<double>>& sample_cov);


double compare_matrices(const std::vector<double> x_paras, 
  const std::vector<std::vector<double>>& pop_cov_mat, 
  std::vector<double>& deviation, std::vector<std::vector<double>>& propose_pop_cov, 
  Eigen::VectorXd& x, Eigen::SparseMatrix<double>& A_mat, const std::vector<int>& index_order);

void optimize_pseudograd_adamw_ver2(std::vector<double>& x0, 
  const std::vector<std::vector<double>>& data_cov_mat, 
  const double lr, const double beta1, const double beta2, const double weight_decay, 
  const int max_step, const double min_error, const double min_abs_error, std::ofstream& log);

void output_matrix(const std::string filename, const std::vector<std::vector<double>>& mat);
void read_mig_matrix(const std::string file_name, std::vector<std::vector<double>>& data_mig_mat);
void initialize_A_mat(const int pop_num, Eigen::SparseMatrix<double>& A_mat, 
  std::vector<int>& index_order);
double softplus(const double x);


#endif
