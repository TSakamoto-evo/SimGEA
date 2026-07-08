#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <string>
#include <sstream>
#include <cmath>
#include <boost/random.hpp>

int main(int argc, char* argv[]){
  int deme_size, gene_len;
  double mig_rate, qtl_mut_rate, neutral_mut_rate;

  if(argc == 7){
    sscanf(argv[2], "%d", &deme_size);
    sscanf(argv[3], "%lf", &mig_rate);
    sscanf(argv[4], "%lf", &qtl_mut_rate);
    sscanf(argv[5], "%lf", &neutral_mut_rate);
    sscanf(argv[6], "%d", &gene_len);
  }else{
    std::cerr << "error" << std::endl;
    std::exit(1);
  }

  std::random_device seed;
  std::mt19937 mt(seed());
  boost::random::mt19937 mt_boost(seed());

  std::vector<double> ss = {0.001, 0.003, 0.01};
  std::shuffle(ss.begin(), ss.end(), mt);
  int selected_num = static_cast<int>(ss.size());

  int gene_num = 200;
  int max_gen = 400000;

  std::vector<double> env;

  {
    // read envfile
    std::ifstream ifs(argv[1]);
    if(!ifs){
      std::cerr << "Fail to open envfile!" << std::endl;
      std::exit(1);
    }

    std::string line;

    while (getline(ifs, line)){
      env.push_back(std::stod(line));
    }
  }

  int pop_num = static_cast<int>(env.size());
  int size = int(std::sqrt(pop_num));

  std::vector<std::vector<int>> selected;
  std::vector<std::vector<int>> neutral;

  std::uniform_int_distribution<> choose_deme(0, pop_num - 1);
  std::poisson_distribution<> neu_mut(2.0 * deme_size * pop_num * neutral_mut_rate * gene_len * gene_num);

  boost::random::binomial_distribution<> det;

  // initialize
  {
    std::vector<int> tmp(pop_num, deme_size);
    for(int i = 0; i < static_cast<int>(ss.size()); i++){
      selected.push_back(tmp);
    }
  }

  for(int gen = 0; gen < max_gen; gen++){
    // selected site
    for(int i = 0; i < selected_num; i++){
      std::vector<int> new_gen(pop_num);

      for(int j = 0; j < pop_num; j++){
        double freq = selected.at(i).at(j) / 2.0 / deme_size;

        // selection
        freq += ss.at(i) * env.at(j) * freq * (1.0 - freq);

        // migration
        if(j / size != 0){
          freq += mig_rate * (selected.at(i).at(size * (j / size - 1) + j % size) - selected.at(i).at(j)) / 2.0 / deme_size;
        }
        if(j / size != size - 1){
          freq += mig_rate * (selected.at(i).at(size * (j / size + 1) + j % size) - selected.at(i).at(j)) / 2.0 / deme_size;
        }
        if(j % size != 0){
          freq += mig_rate * (selected.at(i).at(size * (j / size) + j % size - 1) - selected.at(i).at(j)) / 2.0 / deme_size;
        }
        if(j % size != size - 1){
          freq += mig_rate * (selected.at(i).at(size * (j / size) + j % size + 1) - selected.at(i).at(j)) / 2.0 / deme_size;
        }

        // mutation
        freq += qtl_mut_rate * (1.0 - 2.0 * selected.at(i).at(j) / 2.0 / deme_size);

        det.param(boost::random::binomial_distribution<>::param_type(2 * deme_size, freq));
        new_gen.at(j) = det(mt_boost);
      }

      selected.at(i) = new_gen;
    }

    // neutral site
    for(int i = static_cast<int>(neutral.size()) - 1; i >= 0; i--){
      std::vector<int> new_gen(pop_num);
      int sum_num = 0;

      for(int j = 0; j < pop_num; j++){
        double freq = neutral.at(i).at(j) / 2.0 / deme_size;

        // migration
        if(j / size != 0){
          freq += mig_rate * (neutral.at(i).at(size * (j / size - 1) + j % size) - neutral.at(i).at(j)) / 2.0 / deme_size;
        }
        if(j / size != size - 1){
          freq += mig_rate * (neutral.at(i).at(size * (j / size + 1) + j % size) - neutral.at(i).at(j)) / 2.0 / deme_size;
        }
        if(j % size != 0){
          freq += mig_rate * (neutral.at(i).at(size * (j / size) + j % size - 1) - neutral.at(i).at(j)) / 2.0 / deme_size;
        }
        if(j % size != size - 1){
          freq += mig_rate * (neutral.at(i).at(size * (j / size) + j % size + 1) - neutral.at(i).at(j)) / 2.0 / deme_size;
        }

        det.param(boost::random::binomial_distribution<>::param_type(2 * deme_size, freq));
        new_gen.at(j) = det(mt_boost);
        sum_num += new_gen.at(j);
      }

      if(sum_num == 0 || sum_num == 2 * deme_size * pop_num){
        neutral.erase(neutral.begin() + i);
      }else{
        neutral.at(i) = new_gen;
      }
    }

    // mutation
    int mut_num = neu_mut(mt);
    for(int i = 0; i < mut_num; i++){
      std::vector<int> new_gen(pop_num);
      new_gen.at(choose_deme(mt)) = 1;

      neutral.push_back(new_gen);
    }

    if(gen % 1000 == 0){
      std::cout << gen << "\t" << neutral.size() << std::endl;
    }
  }

  std::ofstream ofs("mean_freq_effect.txt");
  int pos = gene_len * (gene_num / selected_num / 2) + gene_len / 2;
  for(int i = 0; i < selected_num; i++){
    ofs << 0 << "\t" << pos << "\t" << ss.at(i);

    for(const auto& j: selected.at(i)){
      ofs << "\t" << j / 2.0 / deme_size;
    }
    ofs << std::endl;

    pos += gene_len * (gene_num / selected_num);
  }

  std::shuffle(neutral.begin(), neutral.end(), mt);

  std::uniform_int_distribution<> choose_pos(0, gene_num * gene_len - 1);
  std::vector<int> pos_list;

  for(int i = 0; i < static_cast<int>(neutral.size()); i++){
    pos_list.push_back(choose_pos(mt));
  }

  std::sort(pos_list.begin(), pos_list.end());

  std::ofstream ofs2("mean_freq_neutral.txt");
  for(int i = 0; i < static_cast<int>(neutral.size()); i++){
    ofs2 << 0 << "\t" << pos_list.at(i);

    for(const auto& j: neutral.at(i)){
      ofs2 << "\t" << j / 2.0 / deme_size;
    }
    ofs2 << std::endl;
  }
}