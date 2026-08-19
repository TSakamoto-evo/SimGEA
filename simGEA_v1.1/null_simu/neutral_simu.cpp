#include "neutral_simu.hpp"

namespace neusimu{
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

  std::vector<int> get_sample_num(const std::string filename, 
    int& bin_num, std::vector<int>& ret_bin){

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

      int total_num = 0;
      for(int i = 0; i < static_cast<int>(add.size()); i+=2){
        ret.push_back(add.at(i) + add.at(i + 1));
        total_num += add.at(i) + add.at(i + 1);
      }

      bin_num = total_num / 2;

      ret_bin = std::vector<int>(total_num + 1, -1);
      for(int i = 0; i < bin_num; i++){
        ret_bin.at(i + 1) = i;
        ret_bin.at(total_num - (i + 1)) = i;
      }

      return(ret);
    }

    if(ret.empty()){
      std::cerr << "No line in the snpfile!" << std::endl;
      std::exit(1);
    }

    return(ret);
  }

  void read_bin_file(const std::string filename, 
    const int total_num, int& bin_num, std::vector<int>& ret_bin){

    // read binfile
    std::ifstream ifs(filename);
    if(!ifs){
      std::cerr << "Fail to open bin file!" << std::endl;
      std::exit(1);
    }

    std::string line;
    std::vector<double> minor_bin;

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

    if(minor_bin.size() == 0){
      std::cerr << "No valid bin break values" << std::endl;
      std::exit(1);
    }

    std::vector<double> tidy_bins = {minor_bin.at(0)};
    for(int i = 1; i < static_cast<int>(minor_bin.size()); i++){
      if(std::floor(tidy_bins.back() * total_num) != std::floor(minor_bin.at(i) * total_num)){
        tidy_bins.push_back(minor_bin.at(i));
      }
    }

    if(std::floor(tidy_bins.back() * total_num) == std::floor(0.5 * total_num)){
      tidy_bins.pop_back();
    }

    if(tidy_bins.size() == 0){
      std::cerr << "No valid bin break values" << std::endl;
      std::exit(1);
    }
    
    std::ofstream ofs(filename);
    for(const auto& i: tidy_bins){
      ofs << i << std::endl;
    }
    ofs.close();

    bin_num = static_cast<int>(tidy_bins.size());

    ret_bin = std::vector<int>(total_num + 1, -1);
    for(int i = 1; i <= total_num - 1; i++){
      double freq = 1.0 * i / total_num;
      if(freq > 0.5){
        freq = 1.0 - freq;
      }

      for(int j = 0; j < static_cast<int>(tidy_bins.size()); j++){
        if(freq > tidy_bins.at(j) && (j == static_cast<int>(tidy_bins.size() - 1) || freq <= tidy_bins.at(j + 1))){
          ret_bin.at(i) = j;
          break;
        }
      }
    }
  }

  void run_simulations(std::mt19937& mt, const Options& opt, std::ofstream& log, 
    std::vector<int>& sample_num){

    std::vector<int> ret_bin;

    int bin_num;
    sample_num.clear();    
    sample_num = get_sample_num(opt.baypass_file, bin_num, ret_bin);
    int pop_num = static_cast<int>(sample_num.size());

    int total_sample = 0;
    for(int i = 0; i < pop_num; i++){
      total_sample += sample_num.at(i);
    }

    if(opt.sim_mode == 1){
      read_bin_file(opt.bin_file, total_sample, bin_num, ret_bin);
    }

    std::vector<std::vector<double>> mig_mat;
    read_migration_file(opt.mig_file, mig_mat);
    
    if(sample_num.size() != mig_mat.size()){
      std::cerr << "Inconsistent size between the snp file and the migration file" << std::endl;
      std::exit(1);
    }

    int site_per_bin;
    double bin_mut_density;

    if(opt.sim_mode == 0){
      int ac_per_bin = std::min(2 * static_cast<int>(std::floor(opt.window / 2.0 * total_sample)) + 1, total_sample / 2);
      site_per_bin = std::ceil(1.0 * opt.site / ac_per_bin);
      bin_mut_density = opt.mut_density / ac_per_bin;

      log << ac_per_bin << " different allele counts constitute a window\n";
      log << site_per_bin << " mutations will be taken from each allele count\n" << std::endl;
    }else{
      site_per_bin = opt.site;
      bin_mut_density = opt.mut_density;
    }    

    std::uniform_real_distribution<> uni(0.0, 1.0);

    std::vector<int> ini_nums;
    std::vector<std::vector<int>> ini_ind_set(pop_num, std::vector<int>());
    std::vector<int> ini_pos;
    std::vector<std::vector<int>> ini_descendants;

    for(int i = 0; i < pop_num; i++){
      std::vector<int> tmp(pop_num);
      tmp.at(i) = 1;

      ini_nums.push_back(sample_num.at(i));
      for(int j = 0; j < sample_num.at(i); j++){
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
    std::vector<double> branch_time;

    std::vector<double> coal_list(pop_num), mig_list(pop_num);
    
    // pre-run to determine the mutation rate
    std::vector<double> branch_length(bin_num);

    for(int rep = 0; rep < opt.pre_run_tree; rep++){
      nums = ini_nums;
      ind_set = ini_ind_set;
      pos = ini_pos;
      descendants = ini_descendants;
      total_des = std::vector<int>(pos.size(), 1);
      branch_time = std::vector<double>(pos.size(), 0.0);

      int remain = static_cast<int>(pos.size());
      int event = 0;

      double coalescent_rate = 0.0;
      double migration_rate = 0.0;
      double elapsed_time = 0.0;

      for(int i = 0; i < pop_num; i++){
        double n = static_cast<double>(nums.at(i));
        coal_list.at(i) = n * (n - 1) / 2.0;
        coalescent_rate += coal_list.at(i);
        mig_list.at(i) = n * mig_mat.at(i).at(i) * (-0.5);
        migration_rate += mig_list.at(i);
      }

      while(remain > 1){
        event++;

        if(event % 10000 == 0){
          coalescent_rate = 0.0;
          migration_rate = 0.0;

          for(int i = 0; i < pop_num; i++){
            double n = static_cast<double>(nums.at(i));
            coal_list.at(i) = n * (n - 1) / 2.0;
            coalescent_rate += coal_list.at(i);
            mig_list.at(i) = n * mig_mat.at(i).at(i) * (-0.5);
            migration_rate += mig_list.at(i);
          }
        }

        double total_rate = coalescent_rate + migration_rate;
        double time = -1.0 / total_rate * std::log(1.0 - uni(mt));
        elapsed_time += time;

        double u = uni(mt);
        if(u < coalescent_rate / total_rate){
          // coalescent
          double u2 = uni(mt) * coalescent_rate;

          int pop_index = 0;
          while(pop_index < pop_num - 1 && u2 >= coal_list.at(pop_index)){
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

          {
            int bin = ret_bin.at(total_des.at(ind1_index));
            if(bin >= 0){
              branch_length.at(bin) += elapsed_time - branch_time.at(ind1_index);
            }
          }

          {
            int bin = ret_bin.at(total_des.at(ind2_index));
            if(bin >= 0){
              branch_length.at(bin) += elapsed_time - branch_time.at(ind2_index);
            }
          }

          nums.at(pop_index)--;

          for(int i = 0; i < pop_num; i++){
            descendants.at(ind1_index).at(i) += descendants.at(ind2_index).at(i);
          }
          total_des.at(ind1_index) += total_des.at(ind2_index);
          branch_time.at(ind1_index) = elapsed_time;

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
          int pop_index = -1;

          while(pop_index == -1){
            pop_index = 0;
            double u2 = uni(mt) * migration_rate;

            while(pop_index < pop_num - 1 && u2 >= mig_list.at(pop_index)){
              u2 -= mig_list.at(pop_index);
              pop_index++;
            }

            if(nums.at(pop_index) == 0){
              pop_index = -1;
            }
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

    log << "Average branch length\n";

    double shortest_branch = -1.0;
    for(int i = 0; i < bin_num; i++){
      if(branch_length.at(i) <= 0.0){
        if(opt.sim_mode == 0){
          std::cerr << "No positive branch length observed for AC " << i + 1 << ". Increase --prerun or check the model." << std::endl;
          std::exit(1);
        }else{
          std::cerr << "No positive branch length observed for Bin " << i << ". Increase --prerun or check the model." << std::endl;
          std::exit(1);
        }
      }

      if((shortest_branch < 0.0 || shortest_branch > branch_length.at(i))){
        shortest_branch = branch_length.at(i);
      }

      if(opt.sim_mode == 0){
        log << "ac " << i + 1 << ": " << branch_length.at(i) / opt.pre_run_tree << "\n";
      }else{
        log << "bin " << i << ": " << branch_length.at(i) / opt.pre_run_tree << "\n";
      }
    }

    log << "Shortest: " << shortest_branch / opt.pre_run_tree << "\n" << std::endl;
    double per_tree_mut_rate = 1.0 / shortest_branch * opt.pre_run_tree * bin_mut_density;

    std::vector<double> relative_mut(bin_num);
    for(int i = 0; i < bin_num; i++){
      relative_mut.at(i) = shortest_branch / branch_length.at(i);
    }

    std::vector<int> counts(bin_num);

    // simulation
    std::ofstream ofs(opt.simulated_file);
    std::vector<int> sampled_trees(bin_num);
    int valid_bin = bin_num;

    for(long long int rep = 1; rep <= 1e+15; rep++){
      std::vector<int> sampled_or_not(bin_num);

      if(rep % 10000 == 0){
        int min_count = -1;
        for(int i = 0; i < bin_num; i++){
          if(min_count < 0 || min_count > counts.at(i)){
            min_count = counts.at(i);
          }
        }

        log << rep << " trees\t" << min_count << " mutations" << std::endl;
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
        double n = static_cast<double>(nums.at(i));
        coal_list.at(i) = n * (n - 1) / 2.0;
        coalescent_rate += coal_list.at(i);
        mig_list.at(i) = n * mig_mat.at(i).at(i) * (-0.5);
        migration_rate += mig_list.at(i);
        mutation_rate += n * per_tree_mut_rate;
      }

      while(remain > 1){
        event++;

        if(event % 10000 == 0){
          coalescent_rate = 0.0;
          migration_rate = 0.0;
          mutation_rate = 0.0;

          for(int i = 0; i < pop_num; i++){
            double n = static_cast<double>(nums.at(i));
            coal_list.at(i) = n * (n - 1) / 2.0;
            coalescent_rate += coal_list.at(i);
            mig_list.at(i) = n * mig_mat.at(i).at(i) * (-0.5);
            migration_rate += mig_list.at(i);
            mutation_rate += n * per_tree_mut_rate;
          }
        }

        double total_rate = coalescent_rate + migration_rate + mutation_rate;
        
        double u = uni(mt);
        if(u < coalescent_rate / total_rate){
          // coalescent
          double u2 = uni(mt) * coalescent_rate;

          int pop_index = 0;
          while(pop_index < pop_num - 1 && u2 >= coal_list.at(pop_index)){
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
          int pop_index = -1;

          while(pop_index == -1){
            pop_index = 0;
            double u2 = uni(mt) * migration_rate;

            while(pop_index < pop_num - 1 && u2 >= mig_list.at(pop_index)){
              u2 -= mig_list.at(pop_index);
              pop_index++;
            }

            if(nums.at(pop_index) == 0){
              pop_index = -1;
            }
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
          while(pop_index < pop_num - 1 && u2 >= nums.at(pop_index)){
            u2 -= nums.at(pop_index);
            pop_index++;
          }

          int ind1 = std::floor(nums.at(pop_index) * uni(mt));
          if(ind1 == nums.at(pop_index)){
            ind1--;
          }

          int ind_index = ind_set.at(pop_index).at(ind1);
          int bin = ret_bin.at(total_des.at(ind_index));
          
          if(bin >= 0 && counts.at(bin) < site_per_bin && uni(mt) < relative_mut.at(bin)){
            sampled_or_not.at(bin) = 1;

            if(opt.sim_mode == 0){
              ofs << bin + 1;
            }else{
              ofs << bin;
            }
            
            for(int j = 0; j < pop_num; j++){
              ofs << " " << descendants.at(ind_index).at(j) << " " << sample_num.at(j) - descendants.at(ind_index).at(j);
            }
            ofs << std::endl;

            counts.at(bin)++;
            if(counts.at(bin) == site_per_bin){
              if(opt.sim_mode == 0){
                log << rep << " Finish: ac " << bin + 1 << std::endl;
              }else{
                log << rep << " Finish: bin " << bin << std::endl;
              }
              
              valid_bin--;
            }
          }
        }

        if(valid_bin == 0){
          break;
        }
      }

      for(int i = 0; i < bin_num; i++){
        sampled_trees.at(i) += sampled_or_not.at(i);
      }

      if(valid_bin == 0){
        log << "\n" << site_per_bin << " sites from " << rep << " trees\n\n";

        for(int i = 0; i < bin_num; i++){
          if(opt.sim_mode == 0){
            log << "ac " << i + 1 << ": sampled from " << sampled_trees.at(i) << " different trees" << "\n";
          }else{
            log << "bin " << i << ": sampled from " << sampled_trees.at(i) << " different trees" << "\n";
          }
          
        }
        log << std::endl;
        break;
      }
    }

    if(valid_bin != 0){
      std::cerr << "Simulation ended before all allele counts reached the requested number of variants." << std::endl;
      std::exit(1);
    }

    ofs.close();
  }
}