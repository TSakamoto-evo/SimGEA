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
  std::string baypass_file;
  std::string out_prefix = "out";
  std::string log_file = "log.txt";

  std::string version;

  double maf_filter = 0.0;
  double lr = 1e-3;
  double beta1 = 0.9;
  double beta2 = 0.999;
  double weight_decay = 0.01;
  int max_step = 100000;
  double min_error = 1e-6;
  double min_abs_error = 0.01;
  unsigned int seed = 0;
  bool set_seed = 0;
};

bool parse_int_c(const char* s, int& out);
bool parse_double_c(const char* s, double& out);
bool parse_uint_c(const char* s, unsigned int& out);
std::string print_version();


void print_help(const char* prog);
Options parse_args(int argc, char* argv[]);

#endif
