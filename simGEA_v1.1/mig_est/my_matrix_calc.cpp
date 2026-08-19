#include "my_matrix_calc.hpp"

void read_baypass_file(const std::string file_name, 
  std::vector<int>& total_count, 
  std::vector<std::vector<double>>& data_cov_mat, double maf_filter, 
  std::ofstream& log){

  // read snpfile (BayPass style)
  std::ifstream ifs(file_name);
  if(!ifs){
    std::cerr << "Error: fail to open the snpfile!" << std::endl;
    std::exit(1);
  }

  std::string line;

  int line_no = 0;
  double sum_weight = 0.0;
  int snp_no = 0;
  int pop_num = -1;

  maf_filter = std::max(0.0, maf_filter);

  std::vector<double> freq;

  while (getline(ifs, line)){
    std::istringstream iss(line);
    std::string tmp_list;
    std::vector<int> list;

    while(getline(iss, tmp_list, ' ')){
      list.push_back(std::stoi(tmp_list));
    }

    if(pop_num == -1){
      if(list.size() % 2 == 1){
        std::cerr << "Error: Odd number of columns in the snpfile!" << std::endl;
        std::exit(1);
      }

      pop_num = static_cast<int>(list.size()) / 2;
      data_cov_mat = std::vector<std::vector<double>>(pop_num, std::vector<double>(pop_num, 0.0));
      freq = std::vector<double>(pop_num, 0.0);
      total_count = std::vector<int>(pop_num, 0);

      for(int i = 0; i < pop_num; i++){
        total_count[i] = list[2 * i] + list[2 * i + 1];

        if(total_count[i] == 0){
          std::cerr << "Error: No sample in population " << i << " in the snpfile!" << std::endl;
        std::exit(1);
        }
      }

    }else if(2 * pop_num != static_cast<int>(list.size())){
      std::cerr << "Error: # of populations is different at line " << line_no + 1 << " in the snpfile" << std::endl;
      std::exit(1);
    }

    double mean_freq = 0.0;
    for(int i = 0; i < pop_num; i++){
      if(list[2 * i] + list[2 * i + 1] != total_count[i]){
        std::cerr << "Error: # of samples is different at line " << line_no + 1 << " in the snpfile" << std::endl;
        std::exit(1);
      }

      freq[i] = 1.0 * list[2 * i] / (list[2 * i] + list[2 * i + 1]);
      mean_freq += freq[i];
    }
    mean_freq /= pop_num;

    if(mean_freq > maf_filter && 1.0 - mean_freq > maf_filter){
      for(int i = 0; i < pop_num; i++){
        for(int j = i; j < pop_num; j++){
          data_cov_mat[i][j] += (freq[i] - mean_freq) * (freq[j] - mean_freq);
        }
      }

      sum_weight += mean_freq * (1.0 - mean_freq);
      snp_no++;
    }

    line_no++;
  }

  log << "There are " << snp_no << " SNP sites after applying MAF filter of " << maf_filter << "\n" << std::endl;

  if (snp_no == 0 || sum_weight == 0.0) {
    std::cerr << "Error: no SNPs passed the MAF filter (maf=" << maf_filter << ")" << std::endl;
    std::exit(1);
  }

  for(int i = 0; i < pop_num; i++){
    for(int j = i; j < pop_num; j++){
      data_cov_mat[i][j] /= sum_weight;

      if(i != j){
        data_cov_mat[j][i] = data_cov_mat[i][j];
      }
    }
  }
}

void estimate_pop_cov_mat(const std::vector<int>& total_count, 
  const std::vector<std::vector<double>>& sample_cov_mat, 
  std::vector<std::vector<double>>& pop_cov_mat, std::ofstream& log){

  int pop_num = static_cast<int>(sample_cov_mat.size());
  pop_cov_mat = std::vector<std::vector<double>>(pop_num, std::vector<double>(pop_num));
  
  Eigen::MatrixXd A_mat = Eigen::MatrixXd::Zero(pop_num, pop_num);
  Eigen::VectorXd b_vec(pop_num);

  for(int i = 0; i < pop_num; i++){
    b_vec(i) = sample_cov_mat[i][i] - 1;
    A_mat(i, i) = 1.0 / total_count[i] * (1.0 - 2.0 / pop_num) - 1.0;

    for(int j = 0; j < pop_num; j++){
      A_mat(i, j) += 1.0 / pop_num / pop_num / total_count[j];
    }
  }

  Eigen::FullPivLU<Eigen::MatrixXd> LU(A_mat);
  Eigen::VectorXd x = LU.solve(b_vec);

  double last_term = 0.0;
  for(int i = 0; i < pop_num; i++){
    last_term += x(i) / total_count[i];
  }
  last_term /= (pop_num * pop_num);

  Eigen::MatrixXd Pop_mat = Eigen::MatrixXd::Zero(pop_num, pop_num);

  for(int i = 0; i < pop_num; i++){
    Pop_mat(i, i) = 1.0 - x(i);

    for(int j = i + 1; j < pop_num; j++){
      Pop_mat(i, j) = sample_cov_mat[i][j] + x(i) / total_count[i] / pop_num +
        x(j) / total_count[j] / pop_num - last_term;
      Pop_mat(j, i) = Pop_mat(i, j);
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Pop_mat);

  if (solver.info() == Eigen::Success) {
    Eigen::MatrixXd D = solver.eigenvalues().asDiagonal();
    Eigen::MatrixXd V = solver.eigenvectors();
    int negative = 0;

    for(int i = 0; i < pop_num; i++){
      if(D(i, i) < 0.0){
        if(D(i, i) < -1e-6){
          negative++;
        }
        D(i, i) = 0.0;
      }
    }

    if(negative > 0){
      Pop_mat = V * D * V.transpose();
      std::cerr << "Non-small negative eigenvalues in the original Pop matrix: " << 
        negative << " out of " << pop_num << std::endl;
      log << "Non-small negative eigenvalues in the original Pop matrix: " << 
        negative << " out of " << pop_num << std::endl;
    }
  }else{
    std::cerr << "Error: failed in eigen decomposition in the covariance matrix estimation." << std::endl;
    std::exit(1);
  }

  for(int i = 0; i < pop_num; i++){
    for(int j = i; j < pop_num; j++){
      if(i == j){
        pop_cov_mat[i][i] = Pop_mat(i, i);
      }else{
        pop_cov_mat[i][j] = Pop_mat(i, j);
        pop_cov_mat[j][i] = pop_cov_mat[i][j];
      }
    }
  }

  for(int i = 0; i < pop_num; i++){
    double sum_val = 0.0;

    for(int j = 0; j < pop_num; j++){
      if(i != j){
        sum_val += pop_cov_mat[i][j];
      }
    }

    pop_cov_mat[i][i] = -sum_val;
  }
}

void calculate_equilibrium_moment(const std::vector<std::vector<double>>& mig_mat, 
  Eigen::VectorXd& x, Eigen::SparseMatrix<double>& A_mat, const std::vector<int>& index_order){
  
  int pop_num = static_cast<int>(mig_mat.size());
  int mat_size = pop_num * (pop_num + 1) / 2;

  Eigen::VectorXd b_vec(mat_size);

  double* values = A_mat.valuePtr();
  int idx = 0;

  int index1 = 0;
  for(int i = 0; i < pop_num; i++){
    for(int j = i; j < pop_num; j++){
      if(i == j){
        b_vec(index1) = 1.0;
        values[index_order[idx]] = 1.0 - mig_mat[i][i];
        idx++;

        for(int k = 0; k < i; k++){
          values[index_order[idx]] = -mig_mat[i][k];
          idx++;
        }

        for(int k = i + 1; k < pop_num; k++){
          values[index_order[idx]] = -mig_mat[i][k];
          idx++;
        }
      }else{
        b_vec(index1) = 2.0;
        values[index_order[idx]] = - mig_mat[i][i] - mig_mat[j][j];
        idx++;

        for(int k = 0; k < j; k++){
          if(k != i){
            values[index_order[idx]] = -mig_mat[i][k];
            idx++;
          }
        }
        for(int k = j; k < pop_num; k++){
          if(k != i){
            values[index_order[idx]] = -mig_mat[i][k];
            idx++;
          }
        }

        for(int k = 0; k < i; k++){
          if(k != j){
            values[index_order[idx]] = -mig_mat[j][k];
            idx++;
          }
        }
        for(int k = i; k < pop_num; k++){
          if(k != j){
            values[index_order[idx]] = -mig_mat[j][k];
            idx++;
          }
        }
      }

      index1++;
    }
  }

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
  cg.setTolerance(1e-8); 
  cg.compute(A_mat);
  x = cg.solveWithGuess(b_vec, x);

  if (cg.info() != Eigen::Success){
    std::cerr << "CG failed to converge. Use LU factorization." << std::endl;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A_mat);
    solver.factorize(A_mat);
    x = solver.solve(b_vec);
  }
}

void calculate_covariance_structure(const int pop_num, const Eigen::VectorXd& x, 
  std::vector<std::vector<double>>& pop_cov){
  
  pop_cov = std::vector<std::vector<double>>(pop_num, std::vector<double>(pop_num));
  std::vector<double> vec_hi_bar(pop_num);

  int idx = 0;
  for(int i = 0; i < pop_num; i++){
    for(int j = i; j < pop_num; j++){
      if(i == j){
        vec_hi_bar[i] += x(idx);
      }else{
        vec_hi_bar[i] += x(idx);
        vec_hi_bar[j] += x(idx);
      }
      idx++;
    }
  }

  double bar_bar = 0.0;
  for(int i = 0; i < pop_num; i++){
    vec_hi_bar[i] /= pop_num;
    bar_bar += vec_hi_bar[i];
  }
  bar_bar /= pop_num;

  idx = 0;
  for(int i = 0; i < pop_num; i++){
    for(int j = i; j < pop_num; j++){
      if(i == j){
        pop_cov[i][i] = -0.5 * (x(idx) - 2.0 * vec_hi_bar[i] + bar_bar);
        pop_cov[i][i] /= 0.5 * bar_bar;
      }else{
        pop_cov[i][j] = -0.5 * (x(idx) - vec_hi_bar[i] - vec_hi_bar[j] + bar_bar);
        pop_cov[i][j] /= 0.5 * bar_bar;
        pop_cov[j][i] = pop_cov[i][j];
      }
      idx++;
    }
  }

  for(int i = 0; i < pop_num; i++){
    double sum_val = 0.0;
    for(int j = 0; j < pop_num; j++){
      if(i != j){
        sum_val += pop_cov[i][j];
      }
    }
    pop_cov[i][i] = -sum_val;
  }
}

void calculate_covariance_sample_structure(const std::vector<std::vector<double>>& pop_cov, 
  const std::vector<int>& sampled_num, std::vector<std::vector<double>>& sample_cov){

  int pop_num = static_cast<int>(pop_cov.size());
  sample_cov = std::vector<std::vector<double>>(pop_num, std::vector<double>(pop_num));

  if(sampled_num.size() != pop_cov.size()){
    std::cerr << "Error: Length of the sample vector is inconsistent." << std::endl;
    std::exit(1);
  }

  double sum_val = 0.0;
  for(int i = 0; i < pop_num; i++){
    sum_val += (1.0 - pop_cov[i][i]) / sampled_num[i];
  }
  sum_val /= (pop_num * pop_num);

  for(int i = 0; i < pop_num; i++){
    for(int j = i; j < pop_num; j++){
      if(i == j){
        sample_cov[i][i] = 1.0 + (1.0 / sampled_num[i] * (1.0 - 2.0 / pop_num) - 1.0) *
          (1.0 - pop_cov[i][i]) + sum_val;
      }else{
        sample_cov[i][j] = pop_cov[i][j] - 
          1.0 / pop_num * (1.0 - pop_cov[i][i]) / sampled_num[i] - 
          1.0 / pop_num * (1.0 - pop_cov[j][j]) / sampled_num[j] + sum_val;
        sample_cov[j][i] = sample_cov[i][j];
      }
    }
  }
}

double compare_matrices(const std::vector<double> x_paras, 
  const std::vector<std::vector<double>>& pop_cov_mat, 
  std::vector<double>& deviation, std::vector<std::vector<double>>& propose_pop_cov, 
  Eigen::VectorXd& x, Eigen::SparseMatrix<double>& A_mat, const std::vector<int>& index_order){
  
  int pop_num = static_cast<int>(pop_cov_mat.size());
  std::vector<std::vector<double>> propose_mig_mat(pop_num, std::vector<double>(pop_num));
  deviation = std::vector<double>(x_paras.size());

  propose_pop_cov = std::vector<std::vector<double>>(pop_num, std::vector<double>(pop_num));

  int idx = 0;
  for(int i = 0; i < pop_num; i++){
    for(int j = i + 1; j < pop_num; j++){
      propose_mig_mat[i][j] = softplus(x_paras[idx]);
      propose_mig_mat[j][i] = propose_mig_mat[i][j];
      idx++;
    }
  }

  for(int i = 0; i < pop_num; i++){
    double sum_val = 0.0;
    for(int j = 0; j < pop_num; j++){
      if(i != j){
        sum_val += propose_mig_mat[i][j];
      }
    }
    propose_mig_mat[i][i] = -sum_val;
  }

  calculate_equilibrium_moment(propose_mig_mat, x, A_mat, index_order);  
  calculate_covariance_structure(pop_num, x, propose_pop_cov);

  idx = 0;
  for(int i = 0; i < pop_num; i++){
    for(int j = i + 1; j < pop_num; j++){
      deviation[idx] = propose_pop_cov[i][j] - pop_cov_mat[i][j];
      idx++;
    }
  }

  double sum_val = 0.0;
  for(int i = 0; i < pop_num; i++){
    for(int j = 0; j < pop_num; j++){
      if(i != j){
        sum_val += (propose_pop_cov[i][j] - pop_cov_mat[i][j]) * 
          (propose_pop_cov[i][j] - pop_cov_mat[i][j]);
      }
    }
  }

  return(sum_val);
}

void optimize_pseudograd_adamw_ver2(std::vector<double>& x0, 
  const std::vector<std::vector<double>>& data_cov_mat, 
  const double lr, const double beta1, const double beta2, const double weight_decay, 
  const int max_step, const double min_error, const double min_abs_error, std::ofstream& log){

  double eps = 1e-8;
  int vec_len = static_cast<int>(x0.size());
  int pop_num = static_cast<int>(data_cov_mat.size());

  std::vector<double> grad(vec_len);
  std::vector<double> grad_corrected(vec_len);
  std::vector<double> mt(vec_len);
  std::vector<double> vt(vec_len);

  std::vector<double> mt_hat(vec_len);
  std::vector<double> vt_hat(vec_len);

  std::vector<std::vector<double>> pop_cov_mat;

  double previous_loss = 100.0;
  double previous_loss_1000 = 100.0;
  int patience_counter = 0;
  double abs_threshold = 0.0;
  double initial_loss = 0.0;
  bool initialized = 0;

  Eigen::VectorXd x = Eigen::VectorXd::Zero(pop_num * (pop_num + 1) / 2);
  Eigen::SparseMatrix<double> A_mat;
  std::vector<int> index_order;
  initialize_A_mat(pop_num, A_mat, index_order);

  for(int rep = 1; rep <= max_step; rep++){
    double loss = compare_matrices(x0, data_cov_mat, grad, pop_cov_mat, x, A_mat, index_order);

    if(rep == 1){
      initial_loss = loss;
    }

    int idx = 0;
    for(int i = 0; i < pop_num; i++){
      for(int j = i + 1; j < pop_num; j++){
        double sum_grad = (pop_cov_mat[i][i] + pop_cov_mat[j][j] - 2.0 * pop_cov_mat[i][j]) * grad[idx];

        int col_index = i - 1;
        for(int k = 0; k < i; k++){
          if(k != j){
            sum_grad += (pop_cov_mat[j][k] - pop_cov_mat[i][k]) * grad[col_index];
          }
          col_index += (pop_num - k - 2);
        }
        col_index++;
        for(int k = i + 1; k < pop_num; k++){
          if(k != j){
            sum_grad += (pop_cov_mat[j][k] - pop_cov_mat[i][k]) * grad[col_index];
          }
          col_index += 1;
        }

        col_index = j - 1;
        for(int k = 0; k < j; k++){
          if(k != i){
            sum_grad += (pop_cov_mat[i][k] - pop_cov_mat[j][k]) * grad[col_index];
          }
          col_index += (pop_num - k - 2);
        }
        col_index++;
        for(int k = j + 1; k < pop_num; k++){
          if(k != i){
            sum_grad += (pop_cov_mat[i][k] - pop_cov_mat[j][k]) * grad[col_index];
          }
          col_index += 1;
        }

        grad_corrected[idx] = sum_grad / (1.0 + std::exp(-x0[idx]));

        idx++;
      }
    }

    for(int i = 0; i < vec_len; i++){
      mt[i] = beta1 * mt[i] + (1.0 - beta1) * grad_corrected[i];
      vt[i] = beta2 * vt[i] + (1.0 - beta2) * grad_corrected[i] * grad_corrected[i];

      mt_hat[i] = mt[i] / (1.0 - std::pow(beta1, rep));
      vt_hat[i] = vt[i] / (1.0 - std::pow(beta2, rep));

      x0[i] = x0[i] - lr * mt_hat[i] / (std::sqrt(vt_hat[i]) + eps) -
        lr * weight_decay * x0[i];
    }

    if(initialized == 1 && (previous_loss - loss) / (previous_loss + min_error) < min_error &&
      (previous_loss - loss) < abs_threshold){
      patience_counter++;
    }else{
      patience_counter = 0;
    }
    previous_loss = loss;

    if(rep % 1000 == 0){
      log << "Step " << rep << ": " << loss << std::endl;

      if(initialized == 1 && (previous_loss_1000 - loss) / (previous_loss_1000 + min_error) < min_error){
        break;
      }
      previous_loss_1000 = loss;
    }

    if(rep == 1000){
      initialized = 1;
      abs_threshold = (initial_loss - loss) / 1000.0 * min_abs_error;
    }

    if(patience_counter == 10){
      break;
    }
  }
}

void output_matrix(const std::string filename, const std::vector<std::vector<double>>& mat){
  std::ofstream ofs(filename);
  for(int i = 0; i < static_cast<int>(mat.size()); i++){
    ofs << mat[i][0];

    for(int j = 1; j < static_cast<int>(mat[i].size()); j++){
      ofs << " " << mat[i][j];
    }
    ofs << std::endl;
  }
  ofs.close();
}


void read_mig_matrix(const std::string file_name,  
  std::vector<std::vector<double>>& data_mig_mat){

  data_mig_mat.clear();

  // read migfile
  std::ifstream ifs(file_name);
  if(!ifs){
    std::cerr << "Fail to open the migfile!" << std::endl;
    std::exit(1);
  }

  std::string line;

  while (getline(ifs, line)){
    std::istringstream iss(line);
    std::string tmp_list;
    std::vector<double> list;

    while(getline(iss, tmp_list, ' ')){
      list.push_back(std::stod(tmp_list));
    }

    data_mig_mat.push_back(list);
  }
}

void initialize_A_mat(const int pop_num, Eigen::SparseMatrix<double>& A_mat, 
  std::vector<int>& index_order){

  int mat_size = pop_num * (pop_num + 1) / 2;

  std::vector<Eigen::Triplet<double>> triplets;
  double eps = 1e-8;

  std::map<std::pair<int, int>, int> ret_idx;

  int index1 = 0;
  int idx = 0;
  for(int i = 0; i < pop_num; i++){
    for(int j = i; j < pop_num; j++){
      if(i == j){
        triplets.emplace_back(index1, index1, eps);
        ret_idx.emplace(std::make_pair(std::make_pair(index1, index1), idx));
        idx++;

        int col_index = i;
        for(int k = 0; k < i; k++){
          triplets.emplace_back(index1, col_index, eps);
          ret_idx.emplace(std::make_pair(std::make_pair(index1, col_index), idx));
          idx++;
          col_index += (pop_num - k - 1);
        }

        for(int k = i + 1; k < pop_num; k++){
          triplets.emplace_back(index1, index1 - i + k, eps);
          ret_idx.emplace(std::make_pair(std::make_pair(index1, index1 - i + k), idx));
          idx++;
        }
      }else{
        triplets.emplace_back(index1, index1, eps);
        ret_idx.emplace(std::make_pair(std::make_pair(index1, index1), idx));
        idx++;

        int col_index = j;
        for(int k = 0; k < j; k++){
          if(k != i){
            triplets.emplace_back(index1, col_index, eps);
            ret_idx.emplace(std::make_pair(std::make_pair(index1, col_index), idx));
            idx++;
          }
          col_index += (pop_num - k - 1);
        }
        for(int k = j; k < pop_num; k++){
          if(k != i){
            triplets.emplace_back(index1, col_index, eps);
            ret_idx.emplace(std::make_pair(std::make_pair(index1, col_index), idx));
            idx++;
          }
          col_index += 1;
        }

        col_index = i;
        for(int k = 0; k < i; k++){
          if(k != j){
            triplets.emplace_back(index1, col_index, eps);
            ret_idx.emplace(std::make_pair(std::make_pair(index1, col_index), idx));
            idx++;
          }
          col_index += (pop_num - k - 1);
        }
        for(int k = i; k < pop_num; k++){
          if(k != j){
            triplets.emplace_back(index1, col_index, eps);
            ret_idx.emplace(std::make_pair(std::make_pair(index1, col_index), idx));
            idx++;
          }
          col_index += 1;
        }
      }

      index1++;
    }
  }

  if(static_cast<int>(triplets.size()) != idx){
    std::cerr << "Fatal error: failed in constructing the index order vector" << std::endl;
    std::exit(1);
  }

  Eigen::SparseMatrix<double> tmp_A_mat(mat_size, mat_size);
  tmp_A_mat.setFromTriplets(triplets.begin(), triplets.end());

  A_mat = tmp_A_mat;
  A_mat.makeCompressed();

  index_order.clear();
  index_order = std::vector<int>(triplets.size());

  const int* outerIndexPtr = A_mat.outerIndexPtr();

  for (int k = 0; k < A_mat.outerSize(); ++k){
    int count = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(A_mat, k); it; ++it) {
      int row = it.row();
      int col = it.col();

      index_order[ret_idx[{row, col}]] = outerIndexPtr[col] + count;
      count++;
    }
  }
}

double softplus(const double x){
  if(x > 20.0){
    return(x);
  }else if(x < -20.0){
    return(std::exp(x));
  }else{
    return std::log1p(std::exp(x));
  }
}
