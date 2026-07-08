#ifndef OPTIONS
#define OPTIONS

#include <string>
#include <iostream>
#include <vector>
#include <limits>
#include <cstring>
#include <cerrno>
#include <cstdlib>
#include <cmath>

struct Options{
  std::string mig_file;
  std::string baypass_file = "out_data.txt";

  std::string version;

  int site = -1;
  int pre_run_tree = 10000;
  int sample_num = -1;
  double freq_cut = 0.0;
  double exp_mut_per_tree = 0.25;
};

bool parse_int_c(const char* s, int& out);
bool parse_double_c(const char* s, double& out);
std::string print_version();

void print_help(const char* prog);
Options parse_args(int argc, char* argv[]);

#endif
