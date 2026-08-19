#ifndef NEUSIMU
#define NEUSIMU

#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <unordered_set>
#include "options.hpp"

namespace neusimu{
  void read_migration_file(const std::string filename, 
    std::vector<std::vector<double>>& mig_mat);
  std::vector<int> get_sample_num(const std::string filename, 
    int& bin_num, std::vector<int>& ret_bin);
  void read_bin_file(const std::string filename, 
    const int total_num, int& bin_num, std::vector<int>& ret_bin);
  

  void run_simulations(std::mt19937& mt, const Options& opt, std::ofstream& log, 
    std::vector<int>& sample_num);
};

#endif