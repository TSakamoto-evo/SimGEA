#ifndef GEAFUNC
#define GEAFUNC

#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <numeric>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include "options.hpp"

namespace geafunc{
  std::vector<std::vector<double>> read_env_file(const std::string& filename);
  void calculate_null_corr(const std::string filename,
    const std::vector<std::vector<double>>& env_list, 
    const bool use_continuity, 
    std::vector<std::vector<std::vector<double>>>& z_list, 
    const std::vector<int>& sample_num,
    std::vector<std::vector<double>>& ranks_list, 
    std::vector<double>& tie_term_list, 
    std::vector<std::vector<double>>& max_z_list, 
    std::vector<std::vector<double>>& second_max_z_list);
  void prep_tie_groups(const std::vector<double>& env, 
    std::vector<int>& order, std::vector<int>& starts, std::vector<int>& ends);
  double midrank_from_samplenum(const std::vector<int>& order, 
    const std::vector<int>& starts, const std::vector<int>& ends, 
    const std::vector<int>& sample_num, std::vector<double>& ranks);
  void max_u_val(const std::vector<int>& sample_num, 
    const std::vector<double>& ranks, const double tie_term, 
    const bool use_continuity, std::vector<double>& max_z, std::vector<double>& second_max_z);
  void u_from_ranks(const std::vector<int>& c1, const std::vector<int>& c0, 
    const std::vector<double>& ranks, double& ret_u, int& n1, int& n0);
  double z_from_u(const double u, const int n1, const int n0, 
    const double var_u, const bool use_continuity);

  double percentile(const std::vector<std::vector<double>>& z_list, 
    const int low_bin, const int high_bin, const double q);
  double ret_std(const std::vector<std::vector<double>>& vals, 
    const int low_bin, const int high_bin);

  void fit_genpareto(const std::vector<double>& tail_dev, const double max_z, 
    double& ret_c, double& ret_scale, double& min_c, double& max_c, const int grid);
  double best_c_genpareto(const std::vector<double>& tail_dev, const double max_z, 
    const double c, const double prev_scale);
  double gpd_scale_score(const std::vector<double>& tail_dev, const double c, const double scale);
  double gpd_ll(const std::vector<double>& tail_dev, const double c, const double scale);
  
  double calc_kernel_pval(double z_focal, const std::vector<std::vector<double>>& z_list, 
    const double h, const int low_ac, const int high_ac);
  double gpd_sf(const double dev, const double c, const double scale);

  void run_gea_analysis(const Options& opt, std::ofstream& log, 
    const std::vector<int>& sample_num);
};

#endif