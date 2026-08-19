#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <string>
#include <sstream>
#include <cmath>
#include <chrono>
#include "options.hpp"

void read_migration_file(const std::string filename,
  std::vector<std::vector<double>>& mig_mat);

int main(int argc, char* argv[]){
  Options opt = parse_args(argc, argv);

  auto start_time = std::chrono::high_resolution_clock::now();

  int site = opt.site;
  int sample_num = opt.sample_num;
  double freq_cut = opt.freq_cut;
  double exp_mut_per_tree = opt.exp_mut_per_tree;

  std::vector<std::vector<double>> mig_mat;
  read_migration_file(opt.mig_file, mig_mat);

  std::ofstream ofs(opt.baypass_file);

  std::random_device seed;
  std::mt19937 mt(seed());
  std::uniform_real_distribution<> uni(0.0, 1.0);

  int pop_num = static_cast<int>(mig_mat.size());
  std::vector<int> ini_nums;
  std::vector<std::vector<int>> ini_ind_set(pop_num, std::vector<int>());
  std::vector<int> ini_pos;
  std::vector<std::vector<int>> ini_descendants;

  for(int i = 0; i < pop_num; i++){
    std::vector<int> tmp(pop_num);
    tmp.at(i) = 1;

    ini_nums.push_back(sample_num);
    for(int j = 0; j < sample_num; j++){
      ini_ind_set.at(i).push_back(ini_pos.size());
      ini_pos.push_back(i);
      ini_descendants.push_back(tmp);
    }
  }

  std::vector<std::discrete_distribution<>> mig_dests;
  for(int i = 0; i < pop_num; i++){
    std::vector<double> tmp(mig_mat.at(i));
    tmp.at(i) = 0.0;

    std::discrete_distribution<> add(tmp.begin(), tmp.end());
    mig_dests.push_back(add);

    double sum = 0.0;
    for(int j = 0; j < pop_num; j++){
      sum += tmp.at(j);
    }
    mig_mat.at(i).at(i) = -sum;
  }

  std::vector<int> nums;
  std::vector<std::vector<int>> ind_set;
  std::vector<int> pos;
  std::vector<std::vector<int>> descendants;
  std::vector<int> total_des;

  std::vector<double> coal_list(pop_num), mig_list(pop_num);

  // pre-run to determine the mutation rate
  int pre_run_tree = opt.pre_run_tree;
  double branch_length = 0.0;

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
      coal_list.at(i) = nums.at(i) * (nums.at(i) - 1) / 2.0;
      coalescent_rate += coal_list.at(i);
      mig_list.at(i) = nums.at(i) * mig_mat.at(i).at(i) * (-0.5);
      migration_rate += mig_list.at(i);
    }

    while(remain > 1){
      event++;

      if(event % 10000 == 0){
        coalescent_rate = 0.0;
        migration_rate = 0.0;

        for(int i = 0; i < pop_num; i++){
          coal_list.at(i) = nums.at(i) * (nums.at(i) - 1) / 2.0;
          coalescent_rate += coal_list.at(i);
          mig_list.at(i) = nums.at(i) * mig_mat.at(i).at(i) * (-0.5);
          migration_rate += mig_list.at(i);
        }
      }

      double total_rate = coalescent_rate + migration_rate;
      double time = -1.0 / total_rate * std::log(1.0 - uni(mt));

      for(int i = 0; i < pop_num; i++){
        for(int j = 0; j < nums.at(i); j++){
          int ind_index = ind_set.at(i).at(j);

          if(total_des.at(ind_index) > pop_num * sample_num * freq_cut && total_des.at(ind_index) < pop_num * sample_num * (1.0 - freq_cut)){
            branch_length += time;
          }
        }
      }

      double u = uni(mt);
      if(u < coalescent_rate / total_rate){
        // coalescent
        double u2 = uni(mt) * coalescent_rate;

        int pop_index = 0;
        while(pop_index < pop_num - 1 && u2 > coal_list.at(pop_index)){
          u2 -= coal_list.at(pop_index);
          pop_index++;
        }

        int ind1 = std::floor(nums.at(pop_index) * uni(mt));
        if(ind1 == nums.at(pop_index)){
          ind1--;
        }

        int ind2 = std::floor((nums.at(pop_index) - 1) * uni(mt));
        if(ind2 == nums.at(pop_index) - 1){
          ind2--;
        }

        if(ind2 >= ind1){
          ind2++;
        }

        int ind1_index = ind_set.at(pop_index).at(ind1);
        int ind2_index = ind_set.at(pop_index).at(ind2);

        nums.at(pop_index)--;

        for(int i = 0; i < pop_num; i++){
          descendants.at(ind1_index).at(i) += descendants.at(ind2_index).at(i);
        }
        total_des.at(ind1_index) += total_des.at(ind2_index);

        std::swap(ind_set.at(pop_index).at(ind2), ind_set.at(pop_index).back());
        ind_set.at(pop_index).pop_back();
        remain--;

        coal_list.at(pop_index) = nums.at(pop_index) * (nums.at(pop_index) - 1) / 2.0;
        // n*(n+1)/2 - n*(n-1)/2 = n
        coalescent_rate -= nums.at(pop_index);

        mig_list.at(pop_index) = nums.at(pop_index) * mig_mat.at(pop_index).at(pop_index) * (-0.5);
        migration_rate -= mig_mat.at(pop_index).at(pop_index) * (-0.5);
      }else{
        // migration
        double u2 = uni(mt) * migration_rate;

        int pop_index = 0;
        while(pop_index < pop_num - 1 && u2 > mig_list.at(pop_index)){
          u2 -= mig_list.at(pop_index);
          pop_index++;
        }

        int ind1 = std::floor(nums.at(pop_index) * uni(mt));
        if(ind1 == nums.at(pop_index)){
          ind1--;
        }

        int ind1_index = ind_set.at(pop_index).at(ind1);
        int dest = mig_dests.at(pop_index)(mt);

        nums.at(pop_index)--;
        nums.at(dest)++;

        pos.at(ind1_index) = dest;
        
        std::swap(ind_set.at(pop_index).at(ind1), ind_set.at(pop_index).back());
        ind_set.at(pop_index).pop_back();

        ind_set.at(dest).push_back(ind1_index);

        coal_list.at(pop_index) = nums.at(pop_index) * (nums.at(pop_index) - 1) / 2.0;
        coal_list.at(dest) = nums.at(dest) * (nums.at(dest) - 1) / 2.0;
        // n1*(n1-1) / 2 - n1*(n1+1) / 2 + n2*(n2-1) / 2 - (n2-1)*(n2-2)/2
        coalescent_rate += -nums.at(pop_index) + nums.at(dest) - 1;

        mig_list.at(pop_index) = nums.at(pop_index) * mig_mat.at(pop_index).at(pop_index) * (-0.5);
        migration_rate -= mig_mat.at(pop_index).at(pop_index) * (-0.5);

        mig_list.at(dest) = nums.at(dest) * mig_mat.at(dest).at(dest) * (-0.5);
        migration_rate += mig_mat.at(dest).at(dest) * (-0.5);
      }
    }
  }

  branch_length /= pre_run_tree;
  std::cout << branch_length << std::endl;

  double per_tree_mut_rate = 1.0 / branch_length * exp_mut_per_tree;

  //double per_tree_mut_rate = 0.001;
  int total_sampled = 0;

  // simulation
  for(long long int rep = 0; rep < 1e+15; rep++){
    if(rep % 10000 == 0){
      std::cout << rep << std::endl;
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
      coal_list.at(i) = nums.at(i) * (nums.at(i) - 1) / 2.0;
      coalescent_rate += coal_list.at(i);
      mig_list.at(i) = nums.at(i) * mig_mat.at(i).at(i) * (-0.5);
      migration_rate += mig_list.at(i);
      mutation_rate += nums.at(i) * per_tree_mut_rate;
    }

    while(remain > 1){
      event++;

      if(event % 10000 == 0){
        coalescent_rate = 0.0;
        migration_rate = 0.0;
        mutation_rate = 0.0;

        for(int i = 0; i < pop_num; i++){
          coal_list.at(i) = nums.at(i) * (nums.at(i) - 1) / 2.0;
          coalescent_rate += coal_list.at(i);
          mig_list.at(i) = nums.at(i) * mig_mat.at(i).at(i) * (-0.5);
          migration_rate += mig_list.at(i);
          mutation_rate += nums.at(i) * per_tree_mut_rate;
        }
      }

      double total_rate = coalescent_rate + migration_rate + mutation_rate;
      //double time = -1.0 / total_rate * std::log(1.0 - uni(mt));
      
      double u = uni(mt);
      if(u < coalescent_rate / total_rate){
        // coalescent
        double u2 = uni(mt) * coalescent_rate;

        int pop_index = 0;
        while(pop_index < pop_num - 1 && u2 > coal_list.at(pop_index)){
          u2 -= coal_list.at(pop_index);
          pop_index++;
        }

        int ind1 = std::floor(nums.at(pop_index) * uni(mt));
        if(ind1 == nums.at(pop_index)){
          ind1--;
        }

        int ind2 = std::floor((nums.at(pop_index) - 1) * uni(mt));
        if(ind2 == nums.at(pop_index) - 1){
          ind2--;
        }

        if(ind2 >= ind1){
          ind2++;
        }

        int ind1_index = ind_set.at(pop_index).at(ind1);
        int ind2_index = ind_set.at(pop_index).at(ind2);

        nums.at(pop_index)--;

        for(int i = 0; i < pop_num; i++){
          descendants.at(ind1_index).at(i) += descendants.at(ind2_index).at(i);
        }
        total_des.at(ind1_index) += total_des.at(ind2_index);

        std::swap(ind_set.at(pop_index).at(ind2), ind_set.at(pop_index).back());
        ind_set.at(pop_index).pop_back();
        remain--;

        coal_list.at(pop_index) = nums.at(pop_index) * (nums.at(pop_index) - 1) / 2.0;
        coalescent_rate -= nums.at(pop_index);

        mig_list.at(pop_index) = nums.at(pop_index) * mig_mat.at(pop_index).at(pop_index) * (-0.5);
        migration_rate -= mig_mat.at(pop_index).at(pop_index) * (-0.5);

        mutation_rate -= per_tree_mut_rate;
      }else if(u < (coalescent_rate + migration_rate) / total_rate){
        // migration
        double u2 = uni(mt) * migration_rate;

        int pop_index = 0;
        while(pop_index < pop_num - 1 && u2 > mig_list.at(pop_index)){
          u2 -= mig_list.at(pop_index);
          pop_index++;
        }

        int ind1 = std::floor(nums.at(pop_index) * uni(mt));
        if(ind1 == nums.at(pop_index)){
          ind1--;
        }

        int ind1_index = ind_set.at(pop_index).at(ind1);
        int dest = mig_dests.at(pop_index)(mt);

        nums.at(pop_index)--;
        nums.at(dest)++;

        pos.at(ind1_index) = dest;
        
        std::swap(ind_set.at(pop_index).at(ind1), ind_set.at(pop_index).back());
        ind_set.at(pop_index).pop_back();

        ind_set.at(dest).push_back(ind1_index);

        coal_list.at(pop_index) = nums.at(pop_index) * (nums.at(pop_index) - 1) / 2.0;
        coal_list.at(dest) = nums.at(dest) * (nums.at(dest) - 1) / 2.0;
        coalescent_rate += -nums.at(pop_index) + nums.at(dest) - 1;

        mig_list.at(pop_index) = nums.at(pop_index) * mig_mat.at(pop_index).at(pop_index) * (-0.5);
        migration_rate -= mig_mat.at(pop_index).at(pop_index) * (-0.5);

        mig_list.at(dest) = nums.at(dest) * mig_mat.at(dest).at(dest) * (-0.5);
        migration_rate += mig_mat.at(dest).at(dest) * (-0.5);
      }else{
        // mutation
        double u2 = uni(mt) * remain;

        int pop_index = 0;
        while(pop_index < pop_num - 1 && u2 > nums.at(pop_index)){
          u2 -= nums.at(pop_index);
          pop_index++;
        }

        int ind1 = std::floor(nums.at(pop_index) * uni(mt));
        if(ind1 == nums.at(pop_index)){
          ind1--;
        }

        int ind_index = ind_set.at(pop_index).at(ind1);
        int num_des = total_des.at(ind_index);

        if(num_des > pop_num * sample_num * freq_cut && num_des < pop_num * sample_num * (1.0 - freq_cut)){
          ofs << descendants.at(ind_index).at(0) << " " << sample_num - descendants.at(ind_index).at(0);
          for(int j = 1; j < pop_num; j++){
            ofs << " " << descendants.at(ind_index).at(j) << " " << sample_num - descendants.at(ind_index).at(j);
          }
          ofs << std::endl;

          total_sampled++;
          if(total_sampled == site){
            break;
          }
        }
      }
    }

    if(total_sampled == site){
      std::cout << site << " sites from " << rep + 1 << " trees" << std::endl;

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

  std::cout << "Total execution time: " << days << "d " << hours << "h " << minutes << "m " << secs << "s\n";
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
