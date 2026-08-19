#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <filesystem>
#include <sstream>

int main(){
  std::string folder0 = "LE_simu/m=0.04/map0/";
  int sample_num = 0;
  int size = 14;

  std::ofstream ofs(folder0 + "calculate_gst.txt");
  ofs << "data\tidx1\tidx2\tdist\tGst" << std::endl;

  std::ofstream ofs2(folder0 + "site_num.txt");
  ofs2 << "data\tneutral\ttotal" << std::endl;

  for(int data = 1; data <= 50; data++){
    std::string folder = folder0 + "data" + std::to_string(data) + "/";

    std::vector<int> x_pos;
    std::vector<int> y_pos;

    std::string spl = folder + "sampled_pop_list" + std::to_string(sample_num) + ".txt";
    if(std::filesystem::exists(spl)){
      std::ifstream ifs(spl);
      if(!ifs){
        std::cerr << "Fail to open the pos_file!" << std::endl;
        std::exit(1);
      }

      std::string line;

      while (getline(ifs, line)){
        std::istringstream iss(line);
        std::string tmp_list;
        std::vector<std::string> list;

        while(getline(iss, tmp_list, ' ')){
          list.push_back(tmp_list);
        }

        int pop_idx = std::stoi(list.at(0));
        x_pos.push_back((pop_idx - 1) % size);
        y_pos.push_back((pop_idx - 1) / size);
      }
    }else{
      std::cerr << "No pos_file" << std::endl;
    }

    std::vector<bool> neutral;

    std::string ans = folder + "answers" + std::to_string(sample_num) + ".txt";
    if(std::filesystem::exists(ans)){
      std::ifstream ifs(ans);
      if(!ifs){
        std::cerr << "Fail to open the ans_file!" << std::endl;
        std::exit(1);
      }

      std::string line;

      while (getline(ifs, line)){
        std::istringstream iss(line);
        std::string tmp_list;
        std::vector<std::string> list;

        while(getline(iss, tmp_list, ' ')){
          list.push_back(tmp_list);
        }

        if(std::abs(std::stod(list.at(2))) > 1e-8){
          neutral.push_back(0);
        }else{
          neutral.push_back(1);
        }
      }
    }else{
      std::cerr << "No ans_file" << std::endl;
    }

    int pop_num = static_cast<int>(x_pos.size());

    int num_loci = 0;
    int num_line = 0;
    std::vector<std::vector<double>> pi(pop_num, std::vector<double>(pop_num, 0.0));

    std::string fl = folder + "snpsfile" + std::to_string(sample_num) + ".txt";
    if(std::filesystem::exists(fl)){
      std::ifstream ifs(fl);
      if(!ifs){
        std::cerr << "Fail to open the snps_file!" << std::endl;
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

        if(neutral.at(num_line)){
          for(int i = 0; i < pop_num; i++){
            int a0 = list.at(2 * i);
            int a1 = list.at(2 * i + 1);

            double add = 1.0 * a0 * a1 / ((a0 + a1) * (a0 + a1 - 1) / 2);
            pi.at(i).at(i) += add;

            for(int j = i + 1; j < pop_num; j++){
              int b0 = list.at(2 * j);
              int b1 = list.at(2 * j + 1);

              double add = (1.0 * a0 * b1 + 1.0 * a1 * b0) / (1.0 * (a0 + a1) * (b0 + b1));
              pi.at(i).at(j) += add;
            }
          }
          num_loci++;
        }
        
        num_line++;
      }

      ofs2 << data << "\t" << num_loci << "\t" << num_line << std::endl;
    }else{
      std::cerr << "No snps_file" << std::endl;
    }

    for(int i = 0; i < pop_num; i++){
      for(int j = i + 1; j < pop_num; j++){
        double piw = (pi.at(i).at(i) + pi.at(j).at(j)) / 2.0;
        double pit = (pi.at(i).at(i) + pi.at(j).at(j) + 2.0 * pi.at(i).at(j)) / 4.0;

        double gst = (pit - piw) / pit;
        double d = std::sqrt(1.0 * (x_pos.at(i) - x_pos.at(j)) * (x_pos.at(i) - x_pos.at(j)) +
          1.0 * (y_pos.at(i) - y_pos.at(j)) * (y_pos.at(i) - y_pos.at(j)));

        ofs << data << "\t" << i << "\t" << j << "\t" << d << "\t" << gst << std::endl;
      }
    }
  }

  ofs.close();
}