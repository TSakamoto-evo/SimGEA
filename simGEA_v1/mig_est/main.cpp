#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include <chrono>

#include "my_matrix_calc.hpp"
#include "options.hpp"

int main(int argc, char* argv[]){
  Options opt = parse_args(argc, argv);

  std::ofstream log(opt.log_file, std::ios::app);
  log << "MATRIX ESTIMATION\n" << std::endl;

  if(!opt.set_seed){
    std::random_device seed;
    opt.seed = seed();
  }

  // parameter output
  log << "Parameters\n";
  log << "--input " << opt.baypass_file << "\n";
  log << "--out-prefix " << opt.out_prefix << "\n";
  log << "--maf " << opt.maf_filter << "\n";
  log << "--lr " << opt.lr << "\n";
  log << "--beta1 " << opt.beta1 << "\n";
  log << "--beta2 " << opt.beta2 << "\n";
  log << "--weight-decay " << opt.weight_decay << "\n";
  log << "--max-step " << opt.max_step << "\n";
  log << "--min-error " << opt.min_error << "\n";  
  log << "--min-abs-error " << opt.min_abs_error << "\n";  
  log << "--seed " << opt.seed << "\n" << std::endl;
  log << "Version: " << opt.version << "\n" << std::endl;

  auto start_time = std::chrono::high_resolution_clock::now();

  std::vector<std::vector<double>> sample_cov_mat, pop_cov_mat;
  std::vector<int> total_count;

  read_baypass_file(opt.baypass_file, total_count, sample_cov_mat, opt.maf_filter, log);
  estimate_pop_cov_mat(total_count, sample_cov_mat, pop_cov_mat, log);

  output_matrix((opt.out_prefix + "_data_sample_mat.txt").c_str(), sample_cov_mat);
  output_matrix((opt.out_prefix + "_data_pop_mat.txt").c_str(), pop_cov_mat);

  int pop_num = static_cast<int>(total_count.size());

  std::mt19937 mt(opt.seed);
  std::uniform_real_distribution<> initialize(0.0, 1.0);

  std::vector<double> x0(pop_num * (pop_num - 1) / 2);
  for(int i = 0; i < pop_num * (pop_num - 1) / 2; i++){
    x0[i] = initialize(mt);
  }

  log << "start optimization" << std::endl;
  optimize_pseudograd_adamw_ver2(x0, pop_cov_mat, opt.lr, opt.beta1, opt.beta2,
    opt.weight_decay, opt.max_step, opt.min_error, opt.min_abs_error, log);  
  log << "finish optimization\n" << std::endl;

  std::vector<std::vector<double>> final_pop_cov(pop_num, std::vector<double>(pop_num));
  std::vector<std::vector<double>> final_mig_mat(pop_num, std::vector<double>(pop_num));

  int idx = 0;
  for(int i = 0; i < pop_num; i++){
    for(int j = i + 1; j < pop_num; j++){
      final_mig_mat[i][j] = softplus(x0[idx]);
      final_mig_mat[j][i] = final_mig_mat[i][j];
      idx++;
    }
  }

  for(int i = 0; i < pop_num; i++){
    double sum_val = 0.0;
    for(int j = 0; j < pop_num; j++){
      if(i != j){
        sum_val += final_mig_mat[i][j];
      }
    }
    final_mig_mat[i][i] = -sum_val;
  }

  std::vector<std::vector<double>> ret_sample_mat;
  Eigen::VectorXd x = Eigen::VectorXd::Zero(pop_num * (pop_num + 1) / 2);;
  Eigen::SparseMatrix<double> A_mat;
  std::vector<int> index_order;
  initialize_A_mat(pop_num, A_mat, index_order);

  calculate_equilibrium_moment(final_mig_mat, x, A_mat, index_order);
  calculate_covariance_structure(pop_num, x, final_pop_cov);
  calculate_covariance_sample_structure(final_pop_cov, total_count, ret_sample_mat);

  output_matrix((opt.out_prefix + "_est_mig_mat.txt").c_str(), final_mig_mat);
  output_matrix((opt.out_prefix + "_est_cov_mat.txt").c_str(), final_pop_cov);
  output_matrix((opt.out_prefix + "_est_sample_mat.txt").c_str(), ret_sample_mat);

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end_time - start_time;

  int secs = static_cast<int>(elapsed.count());
  int days = secs / (24 * 60 * 60);
  secs -= days * 24 * 60 * 60;
  int hours = secs / (60 * 60);
  secs -= hours * 60 * 60;
  int minutes = secs / 60;
  secs -= minutes * 60;

  log << "Total execution time: " << days << "d " << hours << "h " << minutes << "m " << secs << "s\n";

  log << "END MATRIX ESTIMATION\n\n" << std::endl;

  return 0;
}
