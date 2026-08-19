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
  int mode = 0;
  int sim_mode = 0;

  std::string mig_file;
  std::string baypass_file;
  std::string env_file;
  std::string bin_file;
  std::string simulated_file = "simulated.txt";
  std::string pval_file = "p_vals.txt";
  std::string log_file = "log.txt";

  std::string version;

  int site = -1;
  int pre_run_tree = 10000;
  double mut_density = 0.1;

  double q_tail = 0.995;
  double window = 0.025;

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
