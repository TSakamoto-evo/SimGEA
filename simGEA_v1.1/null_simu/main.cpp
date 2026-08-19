#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <string>
#include <unordered_map>
#include "options.hpp"
#include "neutral_simu.hpp"
#include "gea_func.hpp"

int main(int argc, char* argv[]){
  Options opt = parse_args(argc, argv);

  std::ofstream log(opt.log_file, std::ios::app);
  log << "NEUTRAL SIMULATION\n" << std::endl;

  if(!opt.set_seed){
    std::random_device seed;
    opt.seed = seed();
  }

  // parameter output  
  log << "Parameters\n";
  log << "--mode " << opt.mode << "\n";
  log << "--mig " << opt.mig_file << "\n";
  log << "--data " << opt.baypass_file << "\n";
  log << "--site " << opt.site << "\n";
  log << "--prerun " << opt.pre_run_tree << "\n";
  log << "--mut-density " << opt.mut_density << "\n";
  log << "--output " << opt.simulated_file << "\n\n";

  if(opt.sim_mode == 0){
    log << "--sim-mode " << opt.sim_mode << ": window-based mode\n";
    log << "--win-width " << opt.window << "\n\n";
  }else{
    log << "--sim-mode " << opt.sim_mode << ": bin-based mode\n";
    log << "--bin-file " << opt.bin_file << "\n\n";
  }

  if(opt.mode == 1){
    log << "--env " << opt.env_file << "\n";
    log << "--pval " << opt.pval_file << "\n";
    log << "--q-tail " << opt.q_tail << "\n\n";
  }

  log << "--seed " << opt.seed << "\n";
  log << "Version: " << opt.version << "\n" << std::endl;

  std::mt19937 mt(opt.seed);

  auto start_time = std::chrono::high_resolution_clock::now();

  std::vector<int> sample_num;

  log << "START SIMULATION PART\n" << std::endl;
  neusimu::run_simulations(mt, opt, log, sample_num);
  log << "END SIMULATION PART\n" << std::endl;

  if(opt.mode == 1){
    log << "START GEA PART\n" << std::endl;
    geafunc::run_gea_analysis(opt, log, sample_num);
    log << "END GEA PART\n" << std::endl;
  }

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
  log << "END NEUTRAL SIMULATION\n\n" << std::endl;
}