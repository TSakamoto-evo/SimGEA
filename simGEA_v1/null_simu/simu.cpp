#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <string>
#include <sstream>
#include <cmath>
#include <chrono>
#include <unordered_map>
#include "options.hpp"

void read_migration_file(const std::string filename,
  std::vector<std::vector<double>>& mig_mat);

void read_bin_file(const std::string filename, std::vector<double>& minor_bin);
std::vector<int> get_sample_num(const std::string filename);

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
  log << "--mig " << opt.mig_file << "\n";
  log << "--data " << opt.baypass_file << "\n";
  log << "--bin " << opt.bin_file << "\n";
  log << "--output " << opt.sim_file << "\n";
  log << "--site " << opt.site << "\n";
  log << "--n-pre " << opt.pre_run_tree << "\n";
  log << "--mut-density " << opt.exp_mut_per_tree << "\n";
  log << "--seed " << opt.seed << "\n" << std::endl;
  log << "Version: " << opt.version << "\n" << std::endl;

  auto start_time = std::chrono::high_resolution_clock::now();
  int site = opt.site;
  std::vector<int> sample_num = get_sample_num(opt.baypass_file);

  double exp_mut_per_tree = opt.exp_mut_per_tree;

  std::vector<std::vector<double>> mig_mat;
  read_migration_file(opt.mig_file, mig_mat);

  if(sample_num.size() != mig_mat.size()){
    std::cerr << "Inconsistent size between the snpfile file and the migration file" << std::endl;
    std::exit(1);
  }

  int pop_num = static_cast<int>(mig_mat.size());
  
  int total_sample = 0;
  for(int i = 0; i < pop_num; i++){
    total_sample += sample_num[i];
  }

  std::vector<double> minor_bin;
  read_bin_file(opt.bin_file, minor_bin);

  std::vector<bool> on_off_bin(minor_bin.size() - 1, 1);
  int valid_bin = 0;
  for(int i = 0; i < static_cast<int>(minor_bin.size() - 1); i++){
    if(std::floor(minor_bin[i] * total_sample) == std::floor(minor_bin[i + 1] * total_sample)){
      on_off_bin[i] = 0;
    }else{
      valid_bin++;
    }
  }

  if(valid_bin == 0){
    std::cerr << "No valid bin" << std::endl;
    std::exit(1);
  }

  std::vector<int> ret_bin(total_sample + 1, -1);
  for(int i = 1; i <= total_sample - 1; i++){
    double freq = 1.0 * i / total_sample;
    if(freq > 0.5){
      freq = 1.0 - freq;
    }

    for(int j = 0; j < static_cast<int>(minor_bin.size() - 1); j++){
      if(freq > minor_bin[j] && freq <= minor_bin[j + 1]){
        ret_bin[i] = j;
      }
    }
  }

  std::mt19937 mt(opt.seed);
  std::uniform_real_distribution<> uni(0.0, 1.0);
  
  std::vector<int> ini_nums;
  std::vector<std::vector<int>> ini_ind_set(pop_num, std::vector<int>());
  std::vector<int> ini_pos;
  std::vector<std::vector<int>> ini_descendants;

  for(int i = 0; i < pop_num; i++){
    std::vector<int> tmp(pop_num);
    tmp[i] = 1;

    ini_nums.push_back(sample_num[i]);
    for(int j = 0; j < sample_num[i]; j++){
      ini_ind_set[i].push_back(ini_pos.size());
      ini_pos.push_back(i);
      ini_descendants.push_back(tmp);
    }
  }

  std::ofstream ofs(opt.sim_file);

  std::vector<std::discrete_distribution<>> mig_dests;
  for(int i = 0; i < pop_num; i++){
    std::vector<double> tmp(mig_mat[i]);
    tmp[i] = 0.0;

    std::discrete_distribution<> add(tmp.begin(), tmp.end());
    mig_dests.push_back(add);

    double sum = 0.0;
    for(int j = 0; j < pop_num; j++){
      sum += tmp[j];
    }
    mig_mat[i][i] = -sum;
  }

  std::vector<int> nums;
  std::vector<std::vector<int>> ind_set;
  std::vector<int> pos;
  std::vector<std::vector<int>> descendants;
  std::vector<int> total_des;

  std::vector<double> coal_list(pop_num), mig_list(pop_num);

  // pre-run to determine the mutation rate
  int pre_run_tree = opt.pre_run_tree;
  std::vector<double> branch_length(minor_bin.size() - 1);

  for(int rep = 0; rep < pre_run_tree; rep++){
    nums = ini_nums;
    ind_set = ini_ind_set;
    pos = ini_pos;
    descendants = ini_descendants;
    total_des = std::vector<int>(pos.size(), 1);

    int remain = static_cast<int>(pos.size());
    int event = 0;

    double coalescent_rate = 0.0;
    double migration_rate = 0.0;

    for(int i = 0; i < pop_num; i++){
      coal_list[i] = nums[i] * (nums[i] - 1) / 2.0;
      coalescent_rate += coal_list[i];
      mig_list[i] = nums[i] * mig_mat[i][i] * (-0.5);
      migration_rate += mig_list[i];
    }

    while(remain > 1){
      event++;

      if(event % 10000 == 0){
        coalescent_rate = 0.0;
        migration_rate = 0.0;

        for(int i = 0; i < pop_num; i++){
          coal_list[i] = nums[i] * (nums[i] - 1) / 2.0;
          coalescent_rate += coal_list[i];
          mig_list[i] = nums[i] * mig_mat[i][i] * (-0.5);
          migration_rate += mig_list[i];
        }
      }

      double total_rate = coalescent_rate + migration_rate;
      double time = -1.0 / total_rate * std::log(1.0 - uni(mt));

      for(int i = 0; i < pop_num; i++){
        for(int j = 0; j < nums[i]; j++){
          int ind_index = ind_set[i][j];
          int bin = ret_bin[total_des[ind_index]];

          if(bin >= 0){
            branch_length[bin] += time;
          }
        }
      }

      double u = uni(mt);
      if(u < coalescent_rate / total_rate){
        // coalescent
        double u2 = uni(mt) * coalescent_rate;

        int pop_index = 0;
        while(pop_index < pop_num - 1 && u2 >= coal_list[pop_index]){
          u2 -= coal_list[pop_index];
          pop_index++;
        }

        int ind1 = std::floor(nums[pop_index] * uni(mt));
        if(ind1 == nums[pop_index]){
          ind1--;
        }

        int ind2 = std::floor((nums[pop_index] - 1) * uni(mt));
        if(ind2 == nums[pop_index] - 1){
          ind2--;
        }

        if(ind2 >= ind1){
          ind2++;
        }

        int ind1_index = ind_set[pop_index][ind1];
        int ind2_index = ind_set[pop_index][ind2];

        nums[pop_index]--;

        for(int i = 0; i < pop_num; i++){
          descendants[ind1_index][i] += descendants[ind2_index][i];
        }
        total_des[ind1_index] += total_des[ind2_index];

        std::swap(ind_set[pop_index][ind2], ind_set[pop_index].back());
        ind_set[pop_index].pop_back();
        remain--;

        coal_list[pop_index] = nums[pop_index] * (nums[pop_index] - 1) / 2.0;
        // n*(n+1)/2 - n*(n-1)/2 = n
        coalescent_rate -= nums[pop_index];

        mig_list[pop_index] = nums[pop_index] * mig_mat[pop_index][pop_index] * (-0.5);
        migration_rate -= mig_mat[pop_index][pop_index] * (-0.5);
      }else{
        // migration
        double u2 = uni(mt) * migration_rate;

        int pop_index = 0;
        while(pop_index < pop_num - 1 && u2 >= mig_list[pop_index]){
          u2 -= mig_list[pop_index];
          pop_index++;
        }

        int ind1 = std::floor(nums[pop_index] * uni(mt));
        if(ind1 == nums[pop_index]){
          ind1--;
        }

        int ind1_index = ind_set[pop_index][ind1];
        int dest = mig_dests[pop_index](mt);

        nums[pop_index]--;
        nums[dest]++;

        pos[ind1_index] = dest;
        
        std::swap(ind_set[pop_index][ind1], ind_set[pop_index].back());
        ind_set[pop_index].pop_back();

        ind_set[dest].push_back(ind1_index);

        coal_list[pop_index] = nums[pop_index] * (nums[pop_index] - 1) / 2.0;
        coal_list[dest] = nums[dest] * (nums[dest] - 1) / 2.0;
        // n1*(n1-1) / 2 - n1*(n1+1) / 2 + n2*(n2-1) / 2 - (n2-1)*(n2-2)/2
        coalescent_rate += -nums[pop_index] + nums[dest] - 1;

        mig_list[pop_index] = nums[pop_index] * mig_mat[pop_index][pop_index] * (-0.5);
        migration_rate -= mig_mat[pop_index][pop_index] * (-0.5);

        mig_list[dest] = nums[dest] * mig_mat[dest][dest] * (-0.5);
        migration_rate += mig_mat[dest][dest] * (-0.5);
      }
    }
  }

  log << "\nAverage branch length\n";

  double shortest_branch = -1.0;
  for(int i = 0; i < static_cast<int>(minor_bin.size() - 1); i++){
    if(on_off_bin[i] && (shortest_branch < 0.0 || shortest_branch > branch_length[i])){
      shortest_branch = branch_length[i];
    }

    log << "bin (" << minor_bin[i] << ", " << minor_bin[i+1] << "]: " << branch_length[i] / pre_run_tree << "\n";
  }

  log << "Shortest: " << shortest_branch / pre_run_tree << std::endl;
  double per_tree_mut_rate = 1.0 / shortest_branch * pre_run_tree * exp_mut_per_tree;

  std::vector<double> relative_mut(minor_bin.size() - 1);
  for(int i = 0; i < static_cast<int>(minor_bin.size() - 1); i++){
    if(on_off_bin[i]){
      relative_mut[i] = shortest_branch / branch_length[i];
    }
  }

  std::vector<int> counts(minor_bin.size() - 1);

  // simulation
  for(long long int rep = 1; rep <= 1e+15; rep++){
    if(rep % 10000 == 0){
      int min_count = -1;
      for(int i = 0; i < static_cast<int>(minor_bin.size() - 1); i++){
        if(on_off_bin[i] && (min_count < 0 || min_count > counts[i])){
          min_count = counts[i];
        }
      }

      log << rep << "\t" << min_count << std::endl;
    }

    nums = ini_nums;
    ind_set = ini_ind_set;
    pos = ini_pos;
    descendants = ini_descendants;
    total_des = std::vector<int>(pos.size(), 1);

    int remain = static_cast<int>(pos.size());
    int event = 0;

    double coalescent_rate = 0.0;
    double migration_rate = 0.0;
    double mutation_rate = 0.0;

    for(int i = 0; i < pop_num; i++){
      coal_list[i] = nums[i] * (nums[i] - 1) / 2.0;
      coalescent_rate += coal_list[i];
      mig_list[i] = nums[i] * mig_mat[i][i] * (-0.5);
      migration_rate += mig_list[i];
      mutation_rate += nums[i] * per_tree_mut_rate;
    }

    while(remain > 1){
      event++;

      if(event % 10000 == 0){
        coalescent_rate = 0.0;
        migration_rate = 0.0;
        mutation_rate = 0.0;

        for(int i = 0; i < pop_num; i++){
          coal_list[i] = nums[i] * (nums[i] - 1) / 2.0;
          coalescent_rate += coal_list[i];
          mig_list[i] = nums[i] * mig_mat[i][i] * (-0.5);
          migration_rate += mig_list[i];
          mutation_rate += nums[i] * per_tree_mut_rate;
        }
      }

      double total_rate = coalescent_rate + migration_rate + mutation_rate;
      //double time = -1.0 / total_rate * std::log(1.0 - uni(mt));
      
      double u = uni(mt);
      if(u < coalescent_rate / total_rate){
        // coalescent
        double u2 = uni(mt) * coalescent_rate;

        int pop_index = 0;
        while(pop_index < pop_num - 1 && u2 >= coal_list[pop_index]){
          u2 -= coal_list[pop_index];
          pop_index++;
        }

        int ind1 = std::floor(nums[pop_index] * uni(mt));
        if(ind1 == nums[pop_index]){
          ind1--;
        }

        int ind2 = std::floor((nums[pop_index] - 1) * uni(mt));
        if(ind2 == nums[pop_index] - 1){
          ind2--;
        }

        if(ind2 >= ind1){
          ind2++;
        }

        int ind1_index = ind_set[pop_index][ind1];
        int ind2_index = ind_set[pop_index][ind2];

        nums[pop_index]--;

        for(int i = 0; i < pop_num; i++){
          descendants[ind1_index][i] += descendants[ind2_index][i];
        }
        total_des[ind1_index] += total_des[ind2_index];

        std::swap(ind_set[pop_index][ind2], ind_set[pop_index].back());
        ind_set[pop_index].pop_back();
        remain--;

        coal_list[pop_index] = nums[pop_index] * (nums[pop_index] - 1) / 2.0;
        coalescent_rate -= nums[pop_index];

        mig_list[pop_index] = nums[pop_index] * mig_mat[pop_index][pop_index] * (-0.5);
        migration_rate -= mig_mat[pop_index][pop_index] * (-0.5);

        mutation_rate -= per_tree_mut_rate;
      }else if(u < (coalescent_rate + migration_rate) / total_rate){
        // migration
        double u2 = uni(mt) * migration_rate;

        int pop_index = 0;
        while(pop_index < pop_num - 1 && u2 >= mig_list[pop_index]){
          u2 -= mig_list[pop_index];
          pop_index++;
        }

        int ind1 = std::floor(nums[pop_index] * uni(mt));
        if(ind1 == nums[pop_index]){
          ind1--;
        }

        int ind1_index = ind_set[pop_index][ind1];
        int dest = mig_dests[pop_index](mt);

        nums[pop_index]--;
        nums[dest]++;

        pos[ind1_index] = dest;
        
        std::swap(ind_set[pop_index][ind1], ind_set[pop_index].back());
        ind_set[pop_index].pop_back();

        ind_set[dest].push_back(ind1_index);

        coal_list[pop_index] = nums[pop_index] * (nums[pop_index] - 1) / 2.0;
        coal_list[dest] = nums[dest] * (nums[dest] - 1) / 2.0;
        coalescent_rate += -nums[pop_index] + nums[dest] - 1;

        mig_list[pop_index] = nums[pop_index] * mig_mat[pop_index][pop_index] * (-0.5);
        migration_rate -= mig_mat[pop_index][pop_index] * (-0.5);

        mig_list[dest] = nums[dest] * mig_mat[dest][dest] * (-0.5);
        migration_rate += mig_mat[dest][dest] * (-0.5);
      }else{
        // mutation
        double u2 = uni(mt) * remain;

        int pop_index = 0;
        while(pop_index < pop_num - 1 && u2 >= nums[pop_index]){
          u2 -= nums[pop_index];
          pop_index++;
        }

        int ind1 = std::floor(nums[pop_index] * uni(mt));
        if(ind1 == nums[pop_index]){
          ind1--;
        }

        int ind_index = ind_set[pop_index][ind1];
        int bin = ret_bin[total_des[ind_index]];
        
        if(bin >= 0 && counts[bin] < site && uni(mt) < relative_mut[bin]){
          ofs << bin;
          for(int j = 0; j < pop_num; j++){
            ofs << " " << descendants[ind_index][j] << " " << sample_num[j] - descendants[ind_index][j];
          }
          ofs << std::endl;

          counts[bin]++;
          if(counts[bin] == site){
            log << rep << " Finish: bin " << bin << std::endl;
            valid_bin--;
          }
        }
      }

      if(valid_bin == 0){
        break;
      }
    }

    if(valid_bin == 0){
      log << site << " sites from " << rep << " trees" << std::endl;

      break;
    }
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

void read_migration_file(const std::string filename,
  std::vector<std::vector<double>>& mig_mat){

  // read migfile
  std::ifstream ifs(filename);
  if(!ifs){
    std::cerr << "Fail to open migration file!" << std::endl;
    std::exit(1);
  }

  std::string line;
  int pop_num = -1;

  while (getline(ifs, line)){
    std::istringstream iss(line);
    std::string tmp_list;
    std::vector<double> add;

    while(getline(iss, tmp_list, ' ')){
      add.push_back(std::stod(tmp_list));
    }

    if(pop_num != -1 && pop_num != static_cast<int>(add.size())){
      std::cerr << "The matrix shape is not squared!" << std::endl;
      std::exit(1);
    }
    pop_num = static_cast<int>(add.size());

    mig_mat.push_back(add);
  }

  if(pop_num != static_cast<int>(mig_mat.size())){
    std::cerr << "The matrix shape is not squared!" << std::endl;
    std::exit(1);
  }
}

void read_bin_file(const std::string filename, std::vector<double>& minor_bin){

  // read binfile
  std::ifstream ifs(filename);
  if(!ifs){
    std::cerr << "Fail to open bin file!" << std::endl;
    std::exit(1);
  }

  std::string line;
  minor_bin.clear();

  while (getline(ifs, line)){
    std::istringstream iss(line);

    if(std::stod(line) >= 0.0 && std::stod(line) < 0.5){
      minor_bin.push_back(std::stod(line));
    }else{
      std::cerr << "Some bin values are out of range (ignored)" << std::endl;
    }
  }
  ifs.close();

  std::sort(minor_bin.begin(), minor_bin.end());

  std::ofstream ofs(filename);
  for(const auto& i: minor_bin){
    ofs << i << std::endl;
  }
  ofs.close();

  if(minor_bin.size() == 0){
    std::cerr << "No valid bin value" << std::endl;
    std::exit(1);
  }

  minor_bin.push_back(0.5);
}


std::vector<int> get_sample_num(const std::string filename){
  std::vector<int> ret;

  // read snpfile
  std::ifstream ifs(filename);
  if(!ifs){
    std::cerr << "Fail to open snpfile!" << std::endl;
    std::exit(1);
  }

  std::string line;

  while (getline(ifs, line)){
    std::istringstream iss(line);
    std::string tmp_list;
    std::vector<int> add;

    while(getline(iss, tmp_list, ' ')){
      add.push_back(std::stoi(tmp_list));
    }

    for(int i = 0; i < static_cast<int>(add.size()); i+=2){
      ret.push_back(add[i] + add[i + 1]);
    }

    return(ret);
  }

  std::cerr << "No line in the snpfile!" << std::endl;
  std::exit(1);

  return(ret);
}
