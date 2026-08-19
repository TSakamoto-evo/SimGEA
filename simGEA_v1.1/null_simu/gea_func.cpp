#include "gea_func.hpp"

namespace geafunc{
  std::vector<std::vector<double>> read_env_file(const std::string& filename){
    std::vector<std::vector<double>> ret;

    // read envfile
    std::ifstream ifs(filename);
    if(!ifs){
      std::cerr << "Fail to open env file!" << std::endl;
      std::exit(1);
    }

    std::string line;
    while (getline(ifs, line)){
      std::istringstream iss(line);
      std::string tmp_list;
      std::vector<double> add;

      while(getline(iss, tmp_list, ' ')){
        add.push_back(std::stod(tmp_list));
      }

      ret.push_back(add);
    }

    return(ret);
  }

  void calculate_null_corr(const std::string filename,
    const std::vector<std::vector<double>>& env_list, 
    const bool use_continuity, 
    std::vector<std::vector<std::vector<double>>>& z_list, 
    const std::vector<int>& sample_num,
    std::vector<std::vector<double>>& ranks_list, 
    std::vector<double>& tie_term_list, 
    std::vector<std::vector<double>>& max_z_list, 
    std::vector<std::vector<double>>& second_max_z_list){

    int total_num = 0;
    for(const auto& i: sample_num){
      total_num += i;
    }
    int bin_num = total_num / 2;

    z_list = std::vector<std::vector<std::vector<double>>>(env_list.size(), 
      std::vector<std::vector<double>>(bin_num, std::vector<double>()));
    max_z_list.clear();
    second_max_z_list.clear();

    ranks_list.clear();
    tie_term_list.clear();

    for(size_t i = 0; i < env_list.size(); i++){
      if(sample_num.size() != env_list.at(i).size()){
        std::cerr << "Number of populations is different between the simulated file and envfile!" << std::endl;
        std::exit(1);
      }

      bool identical = 1;
      for(const auto& j: env_list.at(i)){
        if(j != env_list.at(i).at(0)){
          identical = 0;
        }
      }

      if(identical){
        std::cerr << "No variation in environment " << i << std::endl;
        std::exit(1);
      }

      std::vector<int> order, starts, ends;
      prep_tie_groups(env_list.at(i), order, starts, ends);

      std::vector<double> ranks;
      double tie_term = midrank_from_samplenum(order, starts, ends, sample_num, ranks);

      ranks_list.push_back(ranks);
      tie_term_list.push_back(tie_term);

      std::vector<double> max_z, second_max_z;
      max_u_val(sample_num, ranks, tie_term, use_continuity, max_z, second_max_z);
      max_z_list.push_back(max_z);
      second_max_z_list.push_back(second_max_z);
    }

    {
      std::ifstream ifs(filename);
      if(!ifs){
        std::cerr << "Fail to open the simulated file!" << std::endl;
        std::exit(1);
      }

      std::string line;
      while (getline(ifs, line)){
        std::istringstream iss(line);
        std::string tmp_list;
        std::vector<int> list;


        while(getline(iss, tmp_list, ' ')){
          list.push_back(std::stoi(tmp_list));
        }

        int bin_idx = list.at(0) - 1;
        if(bin_idx < 0 || bin_idx >= bin_num){
          std::cerr << "Out of range of bin index in the simulated file" << std::endl;
          std::exit(1);
        }

        std::vector<int> c1, c0;
        int sum_c1 = 0;
        int sum_c0 = 0;

        for(size_t i = 1; i < list.size(); i+=2){
          c1.push_back(list.at(i));
          c0.push_back(list.at(i + 1));
          sum_c1 += list.at(i);
          sum_c0 += list.at(i + 1);

          if(list.at(i) + list.at(i + 1) != sample_num.at(i / 2)){
            std::cerr << "The total number of samples is different" << std::endl;
            std::exit(1);
          }
        }

        if(sum_c0 == 0 || sum_c1 == 0){
          std::cerr << "There is a monomorphic site in the simulated file" << std::endl;
          std::exit(1);
        }

        for(int i = 0; i < static_cast<int>(env_list.size()); i++){
          double u;
          int n1, n0;

          u_from_ranks(c1, c0, ranks_list.at(i), u, n1, n0);
          double var_u = (1.0 * n1 * n0 / 12.0) * (n1 + n0 + 1.0 - tie_term_list.at(i));
          double z_val = z_from_u(u, n1, n0, var_u, use_continuity);

          z_list.at(i).at(bin_idx).push_back(z_val);
          z_list.at(i).at(bin_idx).push_back(-z_val);
        }
      }
    }
  }

  void prep_tie_groups(const std::vector<double>& env, 
    std::vector<int>& order, std::vector<int>& starts, std::vector<int>& ends){

    int n = static_cast<int>(env.size());

    order.resize(n);
    std::iota(order.begin(), order.end(), 0);

    std::sort(order.begin(), order.end(), [&env](int i1, int i2){
      return env.at(i1) < env.at(i2);
    });

    starts.clear();
    ends.clear();
    int start = 0;

    for(int i = 1; i < n; i++){
      if(env.at(order.at(i)) != env.at(order.at(i - 1))){
        starts.push_back(start);
        ends.push_back(i);
        start = i;
      }
    }
    starts.push_back(start);
    ends.push_back(n);
  }

  double midrank_from_samplenum(const std::vector<int>& order, 
    const std::vector<int>& starts, const std::vector<int>& ends, 
    const std::vector<int>& sample_num, std::vector<double>& ranks){
    
    std::vector<int> cumulative_samplenum;
    int cn = 0;
    for(size_t i = 0; i < order.size(); i++){
      cn += sample_num.at(order.at(i));
      cumulative_samplenum.push_back(cn);
    }

    int total_num = 0;
    double tie_term = 0;

    std::vector<double> rank_sorted(order.size());
    for(size_t i = 0; i < starts.size(); i++){
      int left = cumulative_samplenum.at(starts.at(i)) - sample_num.at(order.at(starts.at(i))) + 1;
      int right = cumulative_samplenum.at(ends.at(i) - 1);

      int tie_this_class = 0;
      for(int j = starts.at(i); j < ends.at(i); j++){
        rank_sorted.at(j) = 0.5 * (left + right);
        tie_this_class += sample_num.at(order.at(j));
      }

      total_num += tie_this_class;

      double t = static_cast<double>(tie_this_class);
      tie_term += t * t * t - t;
    }
    tie_term /= 1.0 * total_num * (total_num - 1.0);

    ranks.resize(order.size());
    for(size_t i = 0; i < ranks.size(); i++){
      ranks.at(order.at(i)) = rank_sorted.at(i);
    }

    return(tie_term);
  }

  void max_u_val(const std::vector<int>& sample_num, 
    const std::vector<double>& ranks, const double tie_term, 
    const bool use_continuity, std::vector<double>& max_z, std::vector<double>& second_max_z){
    
    int total_num = 0;
    for(const auto& i: sample_num){
      total_num += i;
    }
    int bin_num = total_num / 2;

    std::vector<int> rank_order(ranks.size());
    std::iota(rank_order.begin(), rank_order.end(), 0);

    std::sort(rank_order.begin(), rank_order.end(), [&ranks](int i1, int i2){
      return ranks.at(i1) < ranks.at(i2);
    });

    max_z = std::vector<double>(bin_num);
    second_max_z = std::vector<double>(bin_num);
    std::vector<double> max_z_gpd(bin_num);

    for(int i = 1; i < total_num; i++){
      int minor_ac = std::min(i, total_num - i);

      std::vector<int> c1(sample_num.size());
      std::vector<int> c0(sample_num);
      int remaining = i;

      for(int j = 0; j < static_cast<int>(sample_num.size()); j++){
        int add = std::min(sample_num.at(rank_order.at(j)), remaining);
        remaining -= add;

        c1.at(rank_order.at(j)) = add;
        c0.at(rank_order.at(j)) -= add;

        if(remaining == 0){
          double ret_u;
          int n1, n0;

          u_from_ranks(c1, c0, ranks, ret_u, n1, n0);
          double var_u = (1.0 * n1 * n0 / 12.0) * (n1 + n0 + 1.0 - tie_term);
          double z = z_from_u(ret_u, n1, n0, var_u, use_continuity);

          bool pattern1 = 0;
          bool pattern2 = 0;

          double z_second, z_third;

          if(add < sample_num.at(rank_order.at(j)) || 
            (j + 1 < static_cast<int>(rank_order.size()) && ranks.at(rank_order.at(j + 1)) == ranks.at(rank_order.at(j)))){
            
            int idx_from = j;
            while(idx_from > 0 && ranks.at(rank_order.at(idx_from)) == ranks.at(rank_order.at(j))){
              idx_from--;
            }

            if(ranks.at(rank_order.at(idx_from)) != ranks.at(rank_order.at(j))){
              std::vector<int> c1_second(c1);
              std::vector<int> c0_second(c0);

              c1_second.at(rank_order.at(idx_from))--;
              c0_second.at(rank_order.at(idx_from))++;

              if(add < sample_num.at(rank_order.at(j))){
                c1_second.at(rank_order.at(j))++;
                c0_second.at(rank_order.at(j))--;
              }else{
                c1_second.at(rank_order.at(j + 1))++;
                c0_second.at(rank_order.at(j + 1))--;
              }

              u_from_ranks(c1_second, c0_second, ranks, ret_u, n1, n0);
              var_u = (1.0 * n1 * n0 / 12.0) * (n1 + n0 + 1.0 - tie_term);
              z_second = z_from_u(ret_u, n1, n0, var_u, use_continuity);

              pattern1 = 1;
            }
          }

          int idx_to = j;
          while(idx_to < static_cast<int>(rank_order.size()) - 1 && ranks.at(rank_order.at(idx_to)) == ranks.at(rank_order.at(j))){
            idx_to++;
          }

          if(ranks.at(rank_order.at(idx_to)) != ranks.at(rank_order.at(j))){
            std::vector<int> c1_second(c1);
            std::vector<int> c0_second(c0);

            c1_second.at(rank_order.at(j))--;
            c0_second.at(rank_order.at(j))++;

            c1_second.at(rank_order.at(idx_to))++;
            c0_second.at(rank_order.at(idx_to))--;

            u_from_ranks(c1_second, c0_second, ranks, ret_u, n1, n0);
            var_u = (1.0 * n1 * n0 / 12.0) * (n1 + n0 + 1.0 - tie_term);
            z_third = z_from_u(ret_u, n1, n0, var_u, use_continuity);

            pattern2 = 1;
          }

          if(pattern1 == 0 && pattern2 == 0){
            z_second = z;
            z_third = z;
          }

          if(pattern1 == 0 || (pattern2 == 1 && std::abs(z_third) > std::abs(z_second))){
            z_second = z_third;
          }

          if(std::abs(z) + (std::abs(z) - std::abs(z_second)) > max_z_gpd.at(minor_ac - 1)){
            max_z_gpd.at(minor_ac - 1) = std::abs(z) + (std::abs(z) - std::abs(z_second));
            max_z.at(minor_ac - 1) = std::abs(z);
            second_max_z.at(minor_ac - 1) = std::abs(z_second);
          }

          break;
        }
      }
    }
  }

  void u_from_ranks(const std::vector<int>& c1, const std::vector<int>& c0, 
    const std::vector<double>& ranks, double& ret_u, int& n1, int& n0){
    
    n0 = 0;
    n1 = 0;
    double r1 = 0.0;

    for(int i = 0; i < static_cast<int>(c1.size()); i++){
      n0 += c0.at(i);
      n1 += c1.at(i);

      r1 += ranks.at(i) * c1.at(i);
    }
    ret_u = r1 - n1 * (n1 + 1.0) / 2.0;
  }

  double z_from_u(const double u, const int n1, const int n0, 
    const double var_u, const bool use_continuity){

    if(var_u <= 0.0){
      return(std::numeric_limits<double>::quiet_NaN());
    }

    double mean_u = 1.0 * n1 * n0 / 2.0;
    double tmp = u - mean_u;

    if(use_continuity){
      double cc = 0.5;
      double tmp2 = std::max(std::abs(tmp) - cc, 0.0);
      
      if(tmp > 0.0){
        tmp = tmp2;
      }else{
        tmp = -tmp2;
      }
    }

    return(tmp / std::sqrt(var_u));
  }


  double percentile(const std::vector<std::vector<double>>& z_list, 
    const int low_bin, const int high_bin, const double q){

    size_t size = 0;
    for(int i = low_bin; i <= high_bin; i++){
      size += z_list.at(i).size();
    }

    std::vector<double> vals;
    vals.reserve(size);

    for(int i = low_bin; i <= high_bin; i++){
      for(const auto& j: z_list.at(i)){
        vals.push_back(j);
      }
    }

    std::sort(vals.begin(), vals.end());

    double index = (vals.size() - 1) * q / 100.0;
    int lower = std::floor(index);
    if(lower == static_cast<int>(vals.size()) - 1){
      lower--;
    }

    double weight = index - lower;

    double ret = vals.at(lower) * (1.0 - weight) + vals.at(lower + 1) * weight;
    return(ret);  
  }

  double ret_std(const std::vector<std::vector<double>>& vals, 
    const int low_bin, const int high_bin){

    size_t size = 0;
    double mean = 0.0;

    for(int i = low_bin; i <= high_bin; i++){
      size += vals.at(i).size();

      for (double x: vals.at(i)){
        mean += x;
      }
    }

    mean /= size;

    double var = 0.0;
    for(int i = low_bin; i <= high_bin; i++){
      for (double x: vals.at(i)){
        double diff = x - mean;
        var += diff * diff;
      }
    }

    var /= static_cast<double>(size - 1);

    return std::sqrt(var);
  }


  void fit_genpareto(const std::vector<double>& tail_dev, const double max_dev, 
    double& ret_c, double& ret_scale, double& min_c, double& max_c, const int grid){

    int best_grid = -1;
    double best_scale = 0.0;
    double best_loglikelihood = -std::numeric_limits<double>::infinity();
    double prev_scale = ret_scale;
    
    for(int i = 0; i <= grid; i++){
      double c = min_c + (max_c - min_c) * i / grid;
      
      prev_scale = best_c_genpareto(tail_dev, max_dev, c, prev_scale);
      double loglikelihood = gpd_ll(tail_dev, c, prev_scale);

      if(loglikelihood > best_loglikelihood){
        best_grid = i;
        best_scale = prev_scale;
        best_loglikelihood = loglikelihood;
      }
    }

    ret_c = min_c + (max_c - min_c) * best_grid / grid;
    ret_scale = best_scale;

    double tmp_min_c = min_c + (max_c - min_c) * std::max(0, best_grid - 1) / grid;
    double tmp_max_c = min_c + (max_c - min_c) * std::min(grid, best_grid + 1) / grid;

    min_c = tmp_min_c;
    max_c = tmp_max_c;
  }

  double best_c_genpareto(const std::vector<double>& tail_dev, const double max_dev, 
    const double shape, const double prev_scale){

    double eps = 1e-8;

    if (shape <= -1.0) {
      std::cerr << "c must be greater than -1" << std::endl;
      std::exit(1);
    }

    if(std::abs(shape) < eps){
      double scale = 0.0;
      for(const auto& j: tail_dev){
        scale += j;
      }
      
      return(scale / tail_dev.size());
    }else{
      double min_scale, max_scale;

      // rough start
      if(prev_scale <= 0.0){
        if(shape >= 0){
          min_scale = eps;
        }else{
          min_scale = -max_dev * shape + eps;
        }

        max_scale = 0.0;
        for(const auto& j: tail_dev){
          max_scale += j;
        }
        max_scale /= tail_dev.size();

        max_scale = std::max(2.0 * min_scale, max_scale);
      }else{
        min_scale = 0.8 * prev_scale;
        max_scale = 1.2 * prev_scale;

        if(shape >= 0 && min_scale < eps){
          min_scale = eps;
        }else if(shape < 0 && min_scale < -max_dev * shape + eps){
          min_scale = -max_dev * shape + eps;
        }else if(gpd_scale_score(tail_dev, shape, min_scale) < 0.0){
          min_scale = (shape >= 0) ? eps : -max_dev * shape + eps;
        }
      }

      if(max_scale <= min_scale){
        max_scale = 2.0 * min_scale + eps;
      }

      // check if score is positive at the minimum possible scale
      if(gpd_scale_score(tail_dev, shape, min_scale) < 0.0){
        return(min_scale);
      }

      // increase max_scale until the score takes positive
      while(gpd_scale_score(tail_dev, shape, max_scale) > 0.0){
        max_scale *= 2;
      }

      double dev;

      for(int iter = 0; iter < 10000; iter++){
        double next_val = (min_scale + max_scale) / 2.0;
        dev = gpd_scale_score(tail_dev, shape, next_val);

        if(dev > 0.0){
          if(min_scale == next_val){
            break;
          }
          min_scale = next_val;
        }else{
          if(max_scale == next_val){
            break;
          }
          max_scale = next_val;
        }

        if(max_scale - min_scale < 1e-16){
          break;
        }
      }

      return((min_scale + max_scale) / 2.0);     
    }
  }

  double gpd_scale_score(const std::vector<double>& tail_dev, const double c, const double scale){
    double ret = -static_cast<double>(tail_dev.size());
    for(const auto& j: tail_dev){
      ret += (1.0 + c) * j / (c * j + scale);
    }

    return(ret);
  }

  double gpd_ll(const std::vector<double>& tail_dev, const double c, const double scale){
    double eps = 1e-8;

    double n = static_cast<double>(tail_dev.size());

    if(std::abs(c) < eps){
      double sum = 0.0;
      for(double j: tail_dev){
        sum += j;
      }
      return(-n * std::log(scale) - sum / scale);
    }else{
      double ret = -n * std::log(scale);
      for(const auto& j: tail_dev){
        ret -= (1.0 + 1.0 / c) * std::log1p(c * j / scale);
      }

      return(ret);
    }
  }

  double calc_kernel_pval(double z_focal, const std::vector<std::vector<double>>& z_list, 
    const double h, const int low_ac, const int high_ac){

    double sum = 0.0;
    size_t size = 0;

    for(int ac = low_ac; ac <= high_ac; ac++){
      size += z_list.at(ac).size();

      for(int i = 0; i < static_cast<int>(z_list.at(ac).size()); i++){
        double tmp = (z_focal - z_list.at(ac).at(i)) / h;
        sum += std::erfc(tmp / std::sqrt(2.0));
      }
    }
    

    sum /= size;
    return(sum);
  }

  double gpd_sf(const double dev, const double c, const double scale){
    if(std::abs(c) < 1e-8){
      return(std::exp(-dev / scale));
    }else{
      return(std::pow(1.0 + c * dev / scale, -1.0 / c));
    }
  }

  void run_gea_analysis(const Options& opt, std::ofstream& log, 
    const std::vector<int>& sample_num){

    int total_num = 0;
    for(const auto& i: sample_num){
      total_num += i;
    }
    int ac_num = total_num / 2;
    int ac_per_bin = std::min(2 * static_cast<int>(std::floor(opt.window / 2.0 * total_num)) + 1, total_num / 2);
    int bin_num = ac_num - ac_per_bin + 1;

    std::vector<int> ret_bin(total_num + 1, -1);

    for(int i = 1; i <= total_num - 1; i++){
      int minor_ac = std::min(i, total_num - i);
      int bin_idx = minor_ac - ac_per_bin / 2 - 1;

      if(bin_idx < 0){
        bin_idx = 0;
      }

      if(bin_idx >= bin_num){
        bin_idx = bin_num - 1;
      }

      ret_bin.at(i) = bin_idx;
    }

    std::vector<std::vector<double>> env_list = read_env_file(opt.env_file);

    std::vector<std::vector<std::vector<double>>> z_list;
    std::vector<std::vector<double>> ranks_list; 
    std::vector<double> tie_term_list;
    std::vector<std::vector<double>> max_z_list;
    std::vector<std::vector<double>> second_max_z_list;

    calculate_null_corr(opt.simulated_file, env_list, 0, z_list, sample_num, 
      ranks_list, tie_term_list, max_z_list, second_max_z_list);
    
    std::vector<std::vector<double>> band_width(env_list.size(), std::vector<double>(ac_num));
    std::vector<std::vector<double>> u_thre(env_list.size(), std::vector<double>(ac_num));
    std::vector<std::vector<double>> cf_thre(env_list.size(), std::vector<double>(ac_num));
    std::vector<std::vector<double>> c_para(env_list.size(), std::vector<double>(ac_num));
    std::vector<std::vector<double>> s_para(env_list.size(), std::vector<double>(ac_num));
    std::vector<std::vector<double>> max_z_gpd(env_list.size(), std::vector<double>(ac_num));

    log << "Start calculating smoothed distribution" << std::endl;

    for(int i = 0; i < static_cast<int>(env_list.size()); i++){
      for(int j = 0; j < ac_num; j++){
        // lowest and highest allele counts used for the estimation of density 
        int low_ac = ret_bin.at(j + 1);
        int high_ac = ret_bin.at(j + 1) + ac_per_bin - 1;

        size_t size = 0;
        for(int ac = low_ac; ac <= high_ac; ac++){
          size += z_list.at(i).at(ac).size();
        }

        double q1 = percentile(z_list.at(i), low_ac, high_ac, 25);
        double q3 = percentile(z_list.at(i), low_ac, high_ac, 75);
        double iqr = q3 - q1;

        double sig_est = std::min(ret_std(z_list.at(i), low_ac, high_ac), iqr / 1.34);
        band_width.at(i).at(j) = 0.9 * sig_est * std::pow(size, -0.2);

        u_thre.at(i).at(j) = percentile(z_list.at(i), low_ac, high_ac, 100 * opt.q_tail);

        double sum = 0.0;
        std::vector<double> tail_dev;
        double max_dev = -1.0;

        for(int ac = low_ac; ac <= high_ac; ac++){
          for(int k = 0; k < static_cast<int>(z_list.at(i).at(ac).size()); k++){
            double tmp = (u_thre.at(i).at(j) - z_list.at(i).at(ac).at(k)) / band_width.at(i).at(j);
            sum += 0.5 * std::erfc(tmp / std::sqrt(2.0));

            if(z_list.at(i).at(ac).at(k) > u_thre.at(i).at(j)){
              tail_dev.push_back(z_list.at(i).at(ac).at(k) - u_thre.at(i).at(j));

              if(tail_dev.back() > max_dev){
                max_dev = tail_dev.back();
              }
            }
          }
        }
        sum /= size;

        if(sum < 0.0){
          sum = 0.0;
        }else if(sum > 1.0){
          sum = 1.0;
        }
        cf_thre.at(i).at(j) = sum;
        max_z_gpd.at(i).at(j) = max_z_list.at(i).at(j) + (max_z_list.at(i).at(j) - second_max_z_list.at(i).at(j)) - u_thre.at(i).at(j);

        if(tail_dev.size() < 250){
          log << "Warning: Tail (q > " << opt.q_tail << ") has only " << tail_dev.size() << " variants for Env " << i << ", AC " << j + 1 << std::endl;
        }

        std::unordered_set<double> unique(tail_dev.begin(), tail_dev.end());

        if(unique.size() > 10){
          double c;
          double scale = 0.0;
          double min_c = -0.999;
          double max_c = 2.0;
          int grid = 100;

          max_dev = std::max(max_dev, max_z_gpd.at(i).at(j));

          fit_genpareto(tail_dev, max_dev, c, scale, min_c, max_c, grid);
          fit_genpareto(tail_dev, max_dev, c, scale, min_c, max_c, grid);
          fit_genpareto(tail_dev, max_dev, c, scale, min_c, max_c, grid);

          c_para.at(i).at(j) = c;
          s_para.at(i).at(j) = scale;
        }else{
          u_thre.at(i).at(j) = u_thre.at(i).at(j) + max_z_gpd.at(i).at(j);
          cf_thre.at(i).at(j) = std::numeric_limits<double>::quiet_NaN();
          c_para.at(i).at(j) = std::numeric_limits<double>::quiet_NaN();
          s_para.at(i).at(j) = std::numeric_limits<double>::quiet_NaN();
        }
      }
    }

    log << "fitted parameters for extreme z_values\n";
    log << "env ac z_thre cf_thre shape scale max_z max_z_gpd" << std::endl;
    for(int i = 0; i < static_cast<int>(env_list.size()); i++){
      for(int j = 0; j < ac_num; j++){
        if(!std::isnan(c_para.at(i).at(j))){
          // gpd was used
          log << i << " " << j + 1 << " " << u_thre.at(i).at(j) << " " << cf_thre.at(i).at(j) << " " << 
            c_para.at(i).at(j) << " " << s_para.at(i).at(j) << " " << max_z_list.at(i).at(j) << " " << u_thre.at(i).at(j) + max_z_gpd.at(i).at(j) << "\n";
        }else{
          log << i << " " << j + 1 << " " << u_thre.at(i).at(j) << " " << cf_thre.at(i).at(j) << " " << 
            c_para.at(i).at(j) << " " << s_para.at(i).at(j) << " " << max_z_list.at(i).at(j) << " " << std::numeric_limits<double>::quiet_NaN() << "\n";
        }
      }
    }

    log << std::endl;

    int num_out_of_range = 0;
    std::ofstream ofs(opt.pval_file);

    std::ifstream ifs(opt.baypass_file);
    if(!ifs){
      std::cerr << "Fail to open data file!" << std::endl;
      std::exit(1);
    }

    std::string line;
    while (getline(ifs, line)){
      std::istringstream iss(line);
      std::string tmp_list;
      std::vector<int> list;

      while(getline(iss, tmp_list, ' ')){
        list.push_back(std::stoi(tmp_list));
      }

      std::vector<int> c1, c0;
      int sum_c1 = 0;

      for(size_t i = 0; i < list.size(); i+=2){
        c1.push_back(list.at(i));
        c0.push_back(list.at(i + 1));
        sum_c1 += list.at(i);

        if(list.at(i) + list.at(i + 1) != sample_num.at(i / 2)){
          std::cerr << "The total number of samples are different" << std::endl;
          std::exit(1);
        }
      }

      int bin_idx = std::min(sum_c1, total_num - sum_c1) - 1;

      if(bin_idx >= 0 && bin_idx < total_num / 2){
        std::vector<double> tmp_pval;
        for(int i = 0; i < static_cast<int>(env_list.size()); i++){
          double u;
          int n1, n0;
          u_from_ranks(c1, c0, ranks_list.at(i), u, n1, n0);
          double var_u = (1.0 * n1 * n0 / 12.0) * (n1 + n0 + 1.0 - tie_term_list.at(i));
          double z_val = z_from_u(u, n1, n0, var_u, 0);

          if(std::abs(z_val) <= u_thre.at(i).at(bin_idx)){
            int low_ac = ret_bin.at(bin_idx + 1);
            int high_ac = ret_bin.at(bin_idx + 1) + ac_per_bin - 1;

            double p_val = calc_kernel_pval(std::abs(z_val), z_list.at(i), band_width.at(i).at(bin_idx), low_ac, high_ac);
            tmp_pval.push_back(p_val);
          }else{
            double p_val = 2.0 * cf_thre.at(i).at(bin_idx) * gpd_sf(std::abs(z_val) - u_thre.at(i).at(bin_idx), 
              c_para.at(i).at(bin_idx), s_para.at(i).at(bin_idx));
            tmp_pval.push_back(p_val);
          }
        }

        ofs << tmp_pval.at(0);
        for(size_t i = 1; i < tmp_pval.size(); i++){
          ofs << " " << tmp_pval.at(i);
        }
        ofs << "\n";
      }else{
        num_out_of_range++;

        ofs << "NA";
        for(int i = 1; i < static_cast<int>(env_list.size()); i++){
          ofs << " NA";
        }
        ofs << "\n";
      }
    }
    ofs.close();

    if(num_out_of_range > 0){
      log << num_out_of_range << " sites are out of null range\n" << std::endl;
    }
  }
}