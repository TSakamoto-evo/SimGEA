functions {
  array [ , ] int construct_idx(int N_patch) {
    array [N_patch ,N_patch] int idx;
    int off_index = N_patch + 1;
    for (i in 1:N_patch) idx[i, i] = i;
    for (i in 1:(N_patch - 1)) {
      for (j in (i + 1):N_patch) {
        idx[i, j] = off_index;
        idx[j, i] = idx[i, j];
        off_index += 1;
      }
    }
    return idx;
  }
  
  matrix Q_rate(matrix M, array [ , ] int idx) {
    int N_patch = rows(M);
    int N_Qrate = (N_patch * (N_patch + 1)) %/% 2;
    int r = N_patch + 1;
    matrix[N_Qrate, N_Qrate] Qrate = rep_matrix(0.0, N_Qrate, N_Qrate); 

    for (i in 1:N_patch) for (j in 1:N_patch) Qrate[i, idx[i, j]] = j==i ? M[i, i] - 1 : M[i, j];
    for (i in 1:(N_patch - 1)) {
      for (j in (i + 1):N_patch) {
        Qrate[r, idx[i, j]] = (M[i, i] + M[j, j]) / 2;
        for (k in 1:N_patch) {
          if(k != j) Qrate[r, i <= k ? idx[i, k] : idx[k, i]] = M[j, k] / 2;
          if(k != i) Qrate[r, j <= k ? idx[j, k] : idx[k, j]] = M[i, k] / 2;
        }
        r += 1;
      }
    }
    return Qrate;
  }
  
  matrix construct_M (vector m_ij, int N_patch){
    matrix[N_patch, N_patch] M = rep_matrix(0.0, N_patch, N_patch);
    int k = 1;
    for(i in 1:(N_patch-1)){
      for(j in (i+1):N_patch){
        M[i,j] = m_ij[k];
        M[j,i] = M[i,j];
        k += 1;
      }
    }
    vector[N_patch] diag_M = M*rep_vector(1,N_patch);
    for(i in 1:N_patch) M[i,i] = -diag_M[i];
    return M;
  }
  
}

data {
  int N_patch;
  vector[(N_patch * (N_patch - 1)) %/% 2] Dij;
}

transformed data{
  int N_Qrate = (N_patch * (N_patch + 1)) %/% 2;
  array [N_patch, N_patch] int idx = construct_idx(N_patch); 
}

parameters{
  vector <lower = 0> [N_Qrate-N_patch] m_ij;
  real <lower = 0> alpha, sigma;
  real gamma0;
}

model{
  vector[N_Qrate - N_patch] mij_tilda = rep_vector(gamma0, N_Qrate - N_patch);
  matrix [N_patch, N_patch] M = construct_M(m_ij, N_patch);
  vector [N_Qrate] t_ij = mdivide_left(Q_rate(M, idx), rep_vector(-1, N_Qrate));
  int pos = 1;
  for(i in 1:(N_patch-1)){
    for(j in (i+1):N_patch){
      Dij[pos] ~ normal(alpha*(4*t_ij[N_patch + pos] - t_ij[i] - t_ij[j]) , sigma);
      pos +=1;
    }
  }
  m_ij ~ exponential(exp(mij_tilda));
  gamma0 ~ normal(0,1);
  [alpha, sigma]'~ exponential(1);
}

generated quantities{
  matrix [N_patch, N_patch] Migration = construct_M(m_ij, N_patch);
}
