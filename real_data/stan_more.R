
### sGMM.onsite ============
sGMM.onsite <- "
data {
  int<lower=0> N;
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K; 
  vector[N] y; 
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;  
  matrix[N, K] B;
  matrix[2, K] B_test;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K] b_spline;
  
  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  //matrix[n_patients, K] r_spline_raw;

  real<lower=0, upper=8> sigma_ref;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  //real<lower=0> sigma_r_spline;
  
  // AR error parameters
  //real<lower=0, upper=0.95> rho_ref;
}

transformed parameters {
  vector[N] mu_y;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  //matrix[n_patients, K] r_spline = sigma_r_spline * r_spline_raw;
  vector[N] mu_fixed = beta0 + beta_bl * x_bl  + beta1 * time + B * b_spline;
  //vector[N] mu_fixed = beta0 + beta1 * time;

  mu_y = mu_fixed + r0[patient_id] + r1[patient_id] .* time;
  //for (i in 1:N) {
  //  mu_y[i] += dot_product(r_spline[patient_id[i], ], B[i, ]); // random spline coef
  //}
}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
 
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,)
  
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  //to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  //sigma_r_spline ~ student_t(1, 0, 5);
  
  // AR error
  //rho_ref ~ uniform(0, 0.95);
  
  // likelihood
  y ~ normal(mu_y, sigma_ref);
  //real sigma_ref_init = sigma_ref / sqrt(1 - rho_ref^2);
  //for (i in 1:N) {
  //  if (ar_idx[i] == 1) {
  //    target += normal_lpdf(y[i] | mu_y[i], sigma_ref_init); //stationary process
  //  } else {
  //    target += normal_lpdf(y[i] | mu_y[i] + rho_ref * (y[i-1] - mu_y[i-1]), sigma_ref);
  //  }
  //}
}

generated quantities {
  //vector[n_patients] mu_onsite_t12 = beta0 + beta1 * 12 + r0 + r1 * 12;
  //array[n_patients] real y_pred_t12;
  //for (ii in 1:n_patients) {
  //  y_pred_t12[ii] = normal_rng(mu_onsite_t12[ii], sigma_ref);
  //}
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  //matrix[n_patients_test, K] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test[2, ], b_spline);
    //real r_spline_contrib = dot_product(r_spline_test[patient_id_test[i], ], B_test[i, ]);
    
    //real mu_fixed_test = beta0 + beta1 * 12;
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12;
    y_pred_t12[i] = mu_y_test[i];
  }
}

"

### sGMM.onsite: v2 with full spline ============
sGMM.onsite.v2 <- "
data {
  int<lower=0> N;
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K1;
  int<lower=0> K2;
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  array[N_obs] int idx_obs;
  array[N_mis] int idx_mis;
  vector[N_obs] y_obs;
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;
  matrix[N, K1] B1;
  matrix[2, K1] B_test1;
  matrix[N, K2] B2;
  matrix[2, K2] B_test2;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K1] b_spline;
  
  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  matrix[n_patients, K2] r_spline_raw;

  real<lower=0, upper=8> sigma_ref;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  real<lower=0> sigma_r_spline;
  vector[N_mis] y_mis;

}

transformed parameters {
  vector[N] y;
  vector[N] mu_y;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  matrix[n_patients, K2] r_spline = sigma_r_spline * r_spline_raw;
  vector[N] mu_fixed = beta0 + beta_bl * x_bl  + beta1 * time + B1 * b_spline;

  y[idx_obs] = y_obs;
  y[idx_mis] = y_mis;
  mu_y = mu_fixed + r0[patient_id] + r1[patient_id] .* time;
  for (i in 1:N) {
    mu_y[i] += dot_product(r_spline[patient_id[i], ], B2[i, ]); // random spline coef
  }
}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
 
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,)
  
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  sigma_r_spline ~ student_t(1, 0, 5);
  
  
  // likelihood
  y ~ normal(mu_y, sigma_ref);
}

generated quantities {
  //vector[n_patients] mu_onsite_t12 = beta0 + beta1 * 12 + r0 + r1 * 12;
  //array[n_patients] real y_pred_t12;
  //for (ii in 1:n_patients) {
  //  y_pred_t12[ii] = normal_rng(mu_onsite_t12[ii], sigma_ref);
  //}
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  vector[K2] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (k in 1:K2) {
    r_spline_test[k] = normal_rng(0, sigma_r_spline);
  }
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test1[2, ], b_spline);
    real r_spline_contrib = dot_product(r_spline_test, B_test2[2, ]);
    
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12 + r_spline_contrib;
    y_pred_t12[i] = mu_y_test[i];
  }
}

"

### sGMM.pool ============
sGMM.pool <- "
data {
  int<lower=0> N; 
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K;             
  int<lower=2> n_cls;  
  int<lower=0> n_ref;
  int<lower=0> n_rem;
  vector[N] y;       
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;
  array[n_ref] int idx_ref;  
  array[n_rem] int idx_rem;  
  matrix[N, K] B;
  matrix[2, K] B_test;
  real<lower=0> var_inflate;
  real<lower=0> a0;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K] b_spline;
  
  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  //matrix[n_patients, K] r_spline_raw;

  real<lower=0, upper=8> sigma_ref;
  real<lower=0, upper=8> sigma_rem;
  //positive_ordered[n_cls-2] sigma_d;
  vector<lower=0, upper=80>[n_cls-2] sigma_d_raw;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  //real<lower=0> sigma_r_spline;
  
  // mixture parameters
  real<lower=0> e0;
  simplex[n_cls] pi; 
}

transformed parameters {
  vector[n_cls-2] sigma_d = sort_asc(sigma_d_raw);
  vector<lower=0>[n_cls] sigma; //remote data
  vector[N] mu_y;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  //matrix[n_patients, K] r_spline = sigma_r_spline * r_spline_raw;
  vector[N] mu_fixed = beta0 + beta_bl * x_bl + beta1 * time + B * b_spline;
  //vector[N] mu_fixed = beta0 + beta1 * time;
  
  sigma[1] = sigma_rem;
  sigma[2] = sqrt(var_inflate) * sigma_ref;
  sigma[3:n_cls] = sqrt(square(sigma[2]) + square(sigma_d));

  mu_y = mu_fixed + r0[patient_id] + r1[patient_id] .* time;
  //for (i in 1:N) {
  //  mu_y[i] += dot_product(r_spline[patient_id[i], ], B[i, ]); // random spline coef
  //}
}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
  sigma_rem ~ student_t(1, 0, 5); // not too large
  sigma_d_raw ~ student_t(3, 0, 5);
  
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,), aka half-Cauchy(5)
  
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  //to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  //sigma_r_spline ~ student_t(1, 0, 5);
  
  
  // mixture weight
  e0 ~ gamma(a0, a0 * n_cls);
  pi ~ dirichlet(rep_vector(e0, n_cls));
  
  //likelihood
  //reference
  y[idx_ref] ~ normal(mu_y[idx_ref], sigma_ref); //ref
  // remote: GMM
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  array[N] int<lower=1, upper=n_cls> z;
  z[idx_ref] = rep_array(1, size(idx_ref));
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    z[i] = categorical_rng(softmax(lps));
  }
  
  //vector[n_patients] mu_onsite_t12 = beta0 + beta1 * 12 + r0 + r1 * 12;
  //array[n_patients] real y_pred_t12;
  //for (ii in 1:n_patients) {
  //  y_pred_t12[ii] = normal_rng(mu_onsite_t12[ii], sigma_ref);
  //}
  
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  //matrix[n_patients_test, K] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test[2, ], b_spline);
    //real r_spline_contrib = dot_product(r_spline_test[patient_id_test[i], ], B_test[i, ]);
    
    //real mu_fixed_test = beta0 + beta1 * 12;
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12;
    y_pred_t12[i] = mu_y_test[i];
  }
  
}

"

### sGMM.pool: v2 with full splines ============
sGMM.pool.v2 <- "
data {
  int<lower=0> N; 
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K1;
  int<lower=0> K2; 
  int<lower=2> n_cls;  
  int<lower=0> n_ref;
  int<lower=0> n_rem;
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  array[N_obs] int idx_obs;
  array[N_mis] int idx_mis;
  vector[N_obs] y_obs;
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;
  array[n_ref] int idx_ref;
  array[n_rem] int idx_rem;
  matrix[N, K1] B1;
  matrix[2, K1] B_test1;
  matrix[N, K2] B2;
  matrix[2, K2] B_test2;
  real<lower=0> var_inflate;
  real<lower=0> a0;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K1] b_spline;

  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  matrix[n_patients, K2] r_spline_raw;

  real<lower=0.8, upper=6> sigma_ref;
  real<lower=0.8, upper=6> sigma_rem;
  vector<lower=0, upper=80>[n_cls-2] sigma_d_raw;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  real<lower=0> sigma_r_spline;

  // mixture parameters
  real<lower=0> e0;
  simplex[n_cls] pi;
  vector[N_mis] y_mis;
}

transformed parameters {
  vector[N] y;
  vector[n_cls-2] sigma_d = sort_asc(sigma_d_raw);
  vector<lower=0>[n_cls] sigma; //remote data
  vector[N] mu_y;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  matrix[n_patients, K2] r_spline = sigma_r_spline * r_spline_raw;
  vector[N] mu_fixed = beta0 + beta_bl * x_bl + beta1 * time + B1 * b_spline;

  y[idx_obs] = y_obs;
  y[idx_mis] = y_mis;
  sigma[1] = sigma_rem;
  sigma[2] = sqrt(var_inflate) * sigma_ref;
  sigma[3:n_cls] = sqrt(square(sigma[2]) + square(sigma_d));

  mu_y = mu_fixed + r0[patient_id] + r1[patient_id] .* time;
  for (i in 1:N) {
    mu_y[i] += dot_product(r_spline[patient_id[i], ], B2[i, ]); // random spline coef
  }
}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
  sigma_rem ~ student_t(1, 0, 5); // not too large
  sigma_d_raw ~ student_t(3, 0, 5);
  
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,), aka half-Cauchy(5)
  
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  sigma_r_spline ~ student_t(1, 0, 5);
  
  
  // mixture weight
  e0 ~ gamma(a0, a0 * n_cls);
  pi ~ dirichlet(rep_vector(e0, n_cls));
  
  //likelihood
  //reference
  y[idx_ref] ~ normal(mu_y[idx_ref], sigma_ref); //ref
  // remote: GMM
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  array[N] int<lower=1, upper=n_cls> z;
  z[idx_ref] = rep_array(1, size(idx_ref));
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    z[i] = categorical_rng(softmax(lps));
  }
  
  
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  vector[K2] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (k in 1:K2) {
    r_spline_test[k] = normal_rng(0, sigma_r_spline);
  }
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test1[2, ], b_spline);
    real r_spline_contrib = dot_product(r_spline_test, B_test2[2, ]);
    
    //real mu_fixed_test = beta0 + beta1 * 12;
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12 + r_spline_contrib;
    y_pred_t12[i] = mu_y_test[i];
  }
  
}

"
### sGMM.pool: v2 with full splines, hetero RE ============
sGMM.pool.v2.hetero <- "
data {
  int<lower=0> N;
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K1;
  int<lower=0> K2;
  int<lower=2> n_cls;  
  int<lower=0> n_ref;
  int<lower=0> n_rem;
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  array[N_obs] int idx_obs;
  array[N_mis] int idx_mis;
  vector[N_obs] y_obs;
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;
  array[n_ref] int idx_ref;
  array[n_rem] int idx_rem;
  matrix[N, K1] B1;
  matrix[2, K1] B_test1;
  matrix[N, K2] B2;
  matrix[2, K2] B_test2;
  real<lower=0> var_inflate;
  real<lower=0> a0;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K1] b_spline;

  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  matrix[n_patients, K2] r_spline_raw;
  vector[n_patients] r0_rem_raw;
  vector[n_patients] r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem_raw;

  real<lower=0.8, upper=6> sigma_ref;
  real<lower=0.8, upper=6> sigma_rem;
  //positive_ordered[n_cls-2] sigma_d;
  vector<lower=0, upper=80>[n_cls-2] sigma_d_raw;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  real<lower=0> sigma_r_spline;
  real<lower=0> sigma_r0_rem;
  real<lower=0> sigma_r1_rem;
  real<lower=0> sigma_r_spline_rem;

  // mixture parameters
  real<lower=0> e0;
  simplex[n_cls] pi;
  vector[N_mis] y_mis;
}

transformed parameters {
  vector[N] y;
  vector[n_cls-2] sigma_d = sort_asc(sigma_d_raw);
  vector<lower=0>[n_cls] sigma; //remote data
  vector[N] mu_y;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  vector[n_patients] r0_rem = sigma_r0_rem * r0_rem_raw;
  vector[n_patients] r1_rem = sigma_r1_rem * r1_rem_raw;
  matrix[n_patients, K2] r_spline = sigma_r_spline * r_spline_raw;
  matrix[n_patients, K2] r_spline_rem  = sigma_r_spline_rem * r_spline_rem_raw;
  vector[N] mu_fixed = beta0 + beta_bl * x_bl + beta1 * time + B1 * b_spline;

  y[idx_obs] = y_obs;
  y[idx_mis] = y_mis;
  sigma[1] = sigma_rem;
  sigma[2] = sqrt(var_inflate) * sigma_ref;
  sigma[3:n_cls] = sqrt(square(sigma[2]) + square(sigma_d));

  for (i in idx_ref) {
    mu_y[i] = mu_fixed[i] + r0[patient_id[i]] + r1[patient_id[i]] * time[i] + dot_product(r_spline[patient_id[i], ], B2[i, ]);
  }
  for (i in idx_rem) {
    mu_y[i] = mu_fixed[i] + r0_rem[patient_id[i]] + r1_rem[patient_id[i]] * time[i] + dot_product(r_spline_rem[patient_id[i], ], B2[i, ]);
  }

}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
  sigma_rem ~ student_t(1, 0, 5); //student_t(1, 0, 5);
  sigma_d_raw ~ student_t(3, 0, 5);
  
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,)
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  sigma_r_spline ~ student_t(1, 0, 5);
  
  r0_rem_raw ~ std_normal();
  r1_rem_raw ~ std_normal();
  to_vector(r_spline_rem_raw) ~ std_normal();
  sigma_r0_rem ~ student_t(1, 0, 5);
  sigma_r1_rem ~ student_t(1, 0, 5);
  sigma_r_spline_rem ~ student_t(1, 0, 5);
  
  // mixture weight
  e0 ~ gamma(a0, a0 * n_cls);
  pi ~ dirichlet(rep_vector(e0, n_cls));

  // likelihood
  //reference
  y[idx_ref] ~ normal(mu_y[idx_ref], sigma_ref);
  
  // remote: GMM
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  array[N] int<lower=1, upper=n_cls> z;

  z[idx_ref] = rep_array(1, size(idx_ref));
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    z[i] = categorical_rng(softmax(lps));
  }
  
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  vector[K2] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (k in 1:K2) {
      r_spline_test[k] = normal_rng(0, sigma_r_spline);
    }
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test1[2, ], b_spline);
    real r_spline_contrib = dot_product(r_spline_test, B_test2[2, ]);
    
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12 + r_spline_contrib;
    y_pred_t12[i] = mu_y_test[i];
  }
  
  
}

"

### sGMM.transfer: v6: ri + delta ============
sGMM.transfer6 <- "
data {
  int<lower=0> N;
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K1;
  int<lower=0> K2;
  int<lower=2> n_cls;  
  int<lower=0> n_ref;
  int<lower=0> n_rem;
  vector[N] y;     
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;  
  array[n_ref] int idx_ref;  
  array[n_rem] int idx_rem;  
  matrix[N, K1] B1;
  matrix[2, K1] B_test1;
  matrix[N, K2] B2;
  matrix[2, K2] B_test2;
  real<lower=0> var_inflate;
  real<lower=0> a0;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K1] b_spline;
  vector[2+K1] dbeta_raw;
  vector<lower=0>[2+K1] lambda_dbeta;
  real<lower=0> tau_dbeta;
  
  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  vector[n_patients] dr0_raw;
  vector[n_patients] dr1_raw;
  matrix[n_patients, K2] r_spline_raw;
  matrix[n_patients, K2] delta_r_spline_raw;
  vector<lower=0>[n_patients] lambda_dr0;
  real<lower=0> tau_dr0;
  vector<lower=0>[n_patients] lambda_dr1;
  real<lower=0> tau_dr1;
  matrix<lower=0>[n_patients, K2] lambda_dr_spline;
  real<lower=0> tau_dr_spline;

  real<lower=0, upper=8> sigma_ref;
  real<lower=0, upper=8> sigma_rem;
  //positive_ordered[n_cls-2] sigma_d;
  vector<lower=0, upper=80>[n_cls-2] sigma_d_raw;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  real<lower=0> sigma_r_spline;
  
  // mixture parameters
  real<lower=0> e0;
  simplex[n_cls] pi; 
}

transformed parameters {
  vector[n_cls-2] sigma_d = sort_asc(sigma_d_raw);
  vector<lower=0>[n_cls] sigma; //remote data
  vector[N] mu_y;
  vector[2+K1] delta_beta = dbeta_raw .* lambda_dbeta * tau_dbeta;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] delta_r0 = dr0_raw .* lambda_dr0 * tau_dr0;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  vector[n_patients] delta_r1 = dr1_raw .* lambda_dr1 * tau_dr1;
  matrix[n_patients, K2] r_spline = sigma_r_spline * r_spline_raw;
  matrix[n_patients, K2] delta_r_spline = delta_r_spline_raw .* lambda_dr_spline * tau_dr_spline;
  vector[N] mu_fixed = beta0 + beta_bl * x_bl + beta1 * time + B1 * b_spline;
  
  sigma[1] = sigma_rem;
  sigma[2] = sqrt(var_inflate) * sigma_ref;
  sigma[3:n_cls] = sqrt(square(sigma[2]) + square(sigma_d));

  mu_y = mu_fixed + r0[patient_id] + r1[patient_id] .* time;
  for (i in 1:N) {
    mu_y[i] += dot_product(r_spline[patient_id[i], ], B2[i, ]); // random spline coef
  }
  mu_y[idx_rem] += delta_beta[1] + delta_beta[2] * time[idx_rem] + B1[idx_rem, ] * delta_beta[3:(2+K1)] + to_vector(delta_r0[patient_id[idx_rem]]) + (to_vector(delta_r1[patient_id[idx_rem]])) .* time[idx_rem] + to_vector(rows_dot_product(delta_r_spline[patient_id[idx_rem], ], B2[idx_rem, ]));
}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
  sigma_rem ~ sigma_prior; //student_t(1, 0, 5);
  sigma_d_raw ~ student_t(3, 0, 5);
  
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,)
  
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  sigma_r_spline ~ student_t(1, 0, 5);
  
  dbeta_raw ~ std_normal();
  lambda_dbeta ~ lambda_prior;
  tau_dbeta ~ tau_prior;
  
  dr0_raw ~ std_normal();
  lambda_dr0 ~ lambda_prior;
  tau_dr0 ~ tau_prior;
  
  dr1_raw ~ std_normal();
  lambda_dr1 ~ lambda_prior;
  tau_dr1 ~ tau_prior;
  
  to_vector(delta_r_spline_raw) ~ std_normal();
  to_vector(lambda_dr_spline) ~ lambda_prior;
  tau_dr_spline ~ tau_prior;
  
  
  // mixture weight
  e0 ~ gamma(a0, a0 * n_cls);
  pi ~ dirichlet(rep_vector(e0, n_cls));

  // likelihood
  //reference
  y[idx_ref] ~ normal(mu_y[idx_ref], sigma_ref);
  
  // remote: GMM
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  array[N] int<lower=1, upper=n_cls> z;

  z[idx_ref] = rep_array(1, size(idx_ref));
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    z[i] = categorical_rng(softmax(lps));
  }
  
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  vector[K2] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (k in 1:K2) {
      r_spline_test[k] = normal_rng(0, sigma_r_spline);
    }
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test1[2, ], b_spline);
    real r_spline_contrib = dot_product(r_spline_test, B_test2[2, ]);
    
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12 + r_spline_contrib;
    y_pred_t12[i] = mu_y_test[i];
  }
  
  
}

"


### sGMM.transfer: v6, seperate ============
sGMM.transfer6.sep <- "
data {
  int<lower=0> N;
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K1;
  int<lower=0> K2;
  int<lower=2> n_cls;  
  int<lower=0> n_ref;
  int<lower=0> n_rem;
  vector[N] y;     
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;  
  array[n_ref] int idx_ref;  
  array[n_rem] int idx_rem;  
  matrix[N, K1] B1;
  matrix[2, K1] B_test1;
  matrix[N, K2] B2;
  matrix[2, K2] B_test2;
  real<lower=0> var_inflate;
  real<lower=0> a0;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K1] b_spline;
  vector[2+K1] delta_beta;
  
  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  matrix[n_patients, K2] r_spline_raw;
  vector[n_patients] r0_rem_raw;
  vector[n_patients] r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem_raw;

  real<lower=0.8, upper=6> sigma_ref;
  real<lower=0.8, upper=6> sigma_rem;
  vector<lower=0, upper=80>[n_cls-2] sigma_d_raw;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  real<lower=0> sigma_r_spline;
  real<lower=0> sigma_r0_rem;
  real<lower=0> sigma_r1_rem;
  real<lower=0> sigma_r_spline_rem;
  
  // mixture parameters
  real<lower=0> e0;
  simplex[n_cls] pi; 
}

transformed parameters {
  vector[n_cls-2] sigma_d = sort_asc(sigma_d_raw);
  vector<lower=0>[n_cls] sigma; //remote data
  vector[N] mu_y;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  matrix[n_patients, K2] r_spline = sigma_r_spline * r_spline_raw;
  vector[n_patients] r0_rem = sigma_r0_rem * r0_rem_raw;
  vector[n_patients] r1_rem = sigma_r1_rem * r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem  = sigma_r_spline_rem * r_spline_rem_raw;
  
  vector[N] mu_fixed = beta0 + beta_bl * x_bl + beta1 * time + B1 * b_spline;
  
  sigma[1] = sigma_rem;
  sigma[2] = sqrt(var_inflate) * sigma_ref;
  sigma[3:n_cls] = sqrt(square(sigma[2]) + square(sigma_d));

  
  for (i in idx_ref) {
    mu_y[i] = mu_fixed[i] + r0[patient_id[i]] + r1[patient_id[i]] .* time[i] + dot_product(r_spline[patient_id[i], ], B2[i, ]);
  }
  for (i in idx_rem) {
    mu_y[i] = mu_fixed[i] + delta_beta[1] + delta_beta[2] * time[i] + B1[i, ] * delta_beta[3:(2+K1)] + r0_rem[patient_id[i]] + r1_rem[patient_id[i]] .* time[i] + dot_product(r_spline_rem[patient_id[i], ], B2[i, ]);
  }
}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
  sigma_rem ~ student_t(1, 0, 5); //student_t(1, 0, 5);
  sigma_d_raw ~ student_t(3, 0, 5);
  
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  delta_beta ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,)
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  sigma_r_spline ~ student_t(1, 0, 5);
  
  r0_rem_raw ~ std_normal();
  r1_rem_raw ~ std_normal();
  to_vector(r_spline_rem_raw) ~ std_normal();
  sigma_r0_rem ~ student_t(1, 0, 5);
  sigma_r1_rem ~ student_t(1, 0, 5);
  sigma_r_spline_rem ~ student_t(1, 0, 5);
  
  // mixture weight
  e0 ~ gamma(a0, a0 * n_cls);
  pi ~ dirichlet(rep_vector(e0, n_cls));

  // likelihood
  //reference
  y[idx_ref] ~ normal(mu_y[idx_ref], sigma_ref);
  
  // remote: GMM
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  array[N] int<lower=1, upper=n_cls> z;

  z[idx_ref] = rep_array(1, size(idx_ref));
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    z[i] = categorical_rng(softmax(lps));
  }
  
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  vector[K2] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (k in 1:K2) {
      r_spline_test[k] = normal_rng(0, sigma_r_spline);
    }
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test1[2, ], b_spline);
    real r_spline_contrib = dot_product(r_spline_test, B_test2[2, ]);
    
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12 + r_spline_contrib;
    y_pred_t12[i] = mu_y_test[i];
  }
  
  
}

"
### sGMM.transfer: v6, semi-seperate ============
sGMM.transfer6.semisep <- "
data {
  int<lower=0> N;
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K1;
  int<lower=0> K2;
  int<lower=2> n_cls;  
  int<lower=0> n_ref;
  int<lower=0> n_rem;
  vector[N] y;     
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;  
  array[n_ref] int idx_ref;  
  array[n_rem] int idx_rem;  
  matrix[N, K1] B1;
  matrix[2, K1] B_test1;
  matrix[N, K2] B2;
  matrix[2, K2] B_test2;
  real<lower=0> var_inflate;
  real<lower=0> a0;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K1] b_spline;
  vector[2+K1] dbeta_raw;
  vector<lower=0>[2+K1] lambda_dbeta;
  real<lower=0> tau_dbeta;
  
  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  matrix[n_patients, K2] r_spline_raw;
  vector[n_patients] r0_rem_raw;
  vector[n_patients] r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem_raw;

  real<lower=0.8, upper=6> sigma_ref;
  real<lower=0.8, upper=6> sigma_rem;
  vector<lower=0, upper=80>[n_cls-2] sigma_d_raw;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  real<lower=0> sigma_r_spline;
  real<lower=0> sigma_r0_rem;
  real<lower=0> sigma_r1_rem;
  real<lower=0> sigma_r_spline_rem;
  
  // mixture parameters
  real<lower=0> e0;
  simplex[n_cls] pi; 
}

transformed parameters {
  vector[n_cls-2] sigma_d = sort_asc(sigma_d_raw);
  vector<lower=0>[n_cls] sigma; //remote data
  vector[N] mu_y;
  vector[2+K1] delta_beta = dbeta_raw .* lambda_dbeta * tau_dbeta;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  matrix[n_patients, K2] r_spline = sigma_r_spline * r_spline_raw;
  vector[n_patients] r0_rem = sigma_r0_rem * r0_rem_raw;
  vector[n_patients] r1_rem = sigma_r1_rem * r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem  = sigma_r_spline_rem * r_spline_rem_raw;
  
  vector[N] mu_fixed = beta0 + beta_bl * x_bl + beta1 * time + B1 * b_spline;
  
  sigma[1] = sigma_rem;
  sigma[2] = sqrt(var_inflate) * sigma_ref;
  sigma[3:n_cls] = sqrt(square(sigma[2]) + square(sigma_d));

  
  for (i in idx_ref) {
    mu_y[i] = mu_fixed[i] + r0[patient_id[i]] + r1[patient_id[i]] .* time[i] + dot_product(r_spline[patient_id[i], ], B2[i, ]);
  }
  for (i in idx_rem) {
    mu_y[i] = mu_fixed[i] + delta_beta[1] + delta_beta[2] * time[i] + B1[i, ] * delta_beta[3:(2+K1)] + r0_rem[patient_id[i]] + r1_rem[patient_id[i]] .* time[i] + dot_product(r_spline_rem[patient_id[i], ], B2[i, ]);
  }
}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
  sigma_rem ~ student_t(1, 0, 5); //student_t(1, 0, 5);
  sigma_d_raw ~ student_t(3, 0, 5);
  
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,)
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  sigma_r_spline ~ student_t(1, 0, 5);
  
  dbeta_raw ~ std_normal();
  lambda_dbeta ~ lambda_prior;
  tau_dbeta ~ tau_prior;
  
  r0_rem_raw ~ std_normal();
  r1_rem_raw ~ std_normal();
  to_vector(r_spline_rem_raw) ~ std_normal();
  sigma_r0_rem ~ student_t(1, 0, 5);
  sigma_r1_rem ~ student_t(1, 0, 5);
  sigma_r_spline_rem ~ student_t(1, 0, 5);
  
  // mixture weight
  e0 ~ gamma(a0, a0 * n_cls);
  pi ~ dirichlet(rep_vector(e0, n_cls));

  // likelihood
  //reference
  y[idx_ref] ~ normal(mu_y[idx_ref], sigma_ref);
  
  // remote: GMM
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  array[N] int<lower=1, upper=n_cls> z;

  z[idx_ref] = rep_array(1, size(idx_ref));
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    z[i] = categorical_rng(softmax(lps));
  }
  
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  vector[K2] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (k in 1:K2) {
      r_spline_test[k] = normal_rng(0, sigma_r_spline);
    }
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test1[2, ], b_spline);
    real r_spline_contrib = dot_product(r_spline_test, B_test2[2, ]);
    
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12 + r_spline_contrib;
    y_pred_t12[i] = mu_y_test[i];
  }
  
  
}

"

### sGMM.transfer: v6: ri, ri*, hs ============
sGMM.transfer6.semisep2 <- "
data {
  int<lower=0> N;
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K1;
  int<lower=0> K2;
  int<lower=2> n_cls;  
  int<lower=0> n_ref;
  int<lower=0> n_rem;
  vector[N] y;     
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;  
  array[n_ref] int idx_ref;  
  array[n_rem] int idx_rem;  
  matrix[N, K1] B1;
  matrix[2, K1] B_test1;
  matrix[N, K2] B2;
  matrix[2, K2] B_test2;
  real<lower=0> var_inflate;
  real<lower=0> a0;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K1] b_spline;
  vector[2+K1] dbeta_raw;
  vector<lower=0>[2+K1] lambda_dbeta;
  real<lower=0> tau_dbeta;
  
  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  matrix[n_patients, K2] r_spline_raw;
  vector[n_patients] r0_rem_raw;
  vector[n_patients] r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem_raw;
  
  vector<lower=0>[n_patients] lambda_r0_rem;
  vector<lower=0>[n_patients] lambda_r1_rem;
  matrix<lower=0>[n_patients, K2] lambda_r_spline_rem;
  real<lower=0> tau_r0_rem;
  real<lower=0> tau_r1_rem;
  real<lower=0> tau_r_spline_rem;
  

  real<lower=0.8, upper=6> sigma_ref;
  real<lower=0.8, upper=6> sigma_rem;
  vector<lower=0, upper=80>[n_cls-2] sigma_d_raw;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  real<lower=0> sigma_r_spline;
  
  // mixture parameters
  real<lower=0> e0;
  simplex[n_cls] pi; 
}

transformed parameters {
  vector[n_cls-2] sigma_d = sort_asc(sigma_d_raw);
  vector<lower=0>[n_cls] sigma; //remote data
  vector[N] mu_y;
  vector[2+K1] delta_beta = dbeta_raw .* lambda_dbeta * tau_dbeta;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  matrix[n_patients, K2] r_spline = sigma_r_spline * r_spline_raw;
  
  vector[n_patients] sigma_r0_rem = sqrt(square(sigma_r0) + square(lambda_r0_rem) * square(tau_r0_rem));
  vector[n_patients] sigma_r1_rem = sqrt(square(sigma_r1) + square(lambda_r1_rem) * square(tau_r1_rem));
  matrix[n_patients, K2] sigma_r_spline_rem = sqrt(square(sigma_r_spline) + square(lambda_r_spline_rem) * square(tau_r_spline_rem));
  
  vector[n_patients] r0_rem = sigma_r0_rem .* r0_rem_raw;
  vector[n_patients] r1_rem = sigma_r1_rem .* r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem  = sigma_r_spline_rem .* r_spline_rem_raw;
  
  vector[N] mu_fixed = beta0 + beta_bl * x_bl + beta1 * time + B1 * b_spline;
  
  sigma[1] = sigma_rem;
  sigma[2] = sqrt(var_inflate) * sigma_ref;
  sigma[3:n_cls] = sqrt(square(sigma[2]) + square(sigma_d));

  
  for (i in idx_ref) {
    mu_y[i] = mu_fixed[i] + r0[patient_id[i]] + r1[patient_id[i]] * time[i] + dot_product(r_spline[patient_id[i], ], B2[i, ]);
  }
  for (i in idx_rem) {
    mu_y[i] = mu_fixed[i] + delta_beta[1] + delta_beta[2] * time[i] + B1[i, ] * delta_beta[3:(2+K1)] + r0_rem[patient_id[i]] + r1_rem[patient_id[i]] * time[i] + dot_product(r_spline_rem[patient_id[i], ], B2[i, ]);
  }
}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
  sigma_rem ~ student_t(1, 0, 5); //student_t(1, 0, 5);
  sigma_d_raw ~ student_t(3, 0, 5);
  
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,)
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  sigma_r_spline ~ student_t(1, 0, 5);
  
  dbeta_raw ~ std_normal();
  lambda_dbeta ~ lambda_prior;
  tau_dbeta ~ tau_prior;
  
  r0_rem_raw ~ std_normal();
  r1_rem_raw ~ std_normal();
  to_vector(r_spline_rem_raw) ~ std_normal();
  
  lambda_r0_rem ~ lambda_prior;
  lambda_r1_rem ~ lambda_prior;
  to_vector(lambda_r_spline_rem) ~ lambda_prior;

  tau_r0_rem ~ tau_prior;
  tau_r1_rem ~ tau_prior;
  tau_r_spline_rem ~ tau_prior;
  
  // mixture weight
  e0 ~ gamma(a0, a0 * n_cls);
  pi ~ dirichlet(rep_vector(e0, n_cls));

  // likelihood
  //reference
  y[idx_ref] ~ normal(mu_y[idx_ref], sigma_ref);
  
  // remote: GMM
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  array[N] int<lower=1, upper=n_cls> z;

  z[idx_ref] = rep_array(1, size(idx_ref));
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    z[i] = categorical_rng(softmax(lps));
  }
  
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  vector[K2] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (k in 1:K2) {
      r_spline_test[k] = normal_rng(0, sigma_r_spline);
    }
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test1[2, ], b_spline);
    real r_spline_contrib = dot_product(r_spline_test, B_test2[2, ]);
    
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12 + r_spline_contrib;
    y_pred_t12[i] = mu_y_test[i];
  }
  
  
}

"

### sGMM.transfer: v6: ri, ri*, group hs ============
sGMM.transfer6.semisep3 <- "
data {
  int<lower=0> N;
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K1;
  int<lower=0> K2;
  int<lower=2> n_cls;  
  int<lower=0> n_ref;
  int<lower=0> n_rem;
  vector[N] y;     
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;  
  array[n_ref] int idx_ref;  
  array[n_rem] int idx_rem;  
  matrix[N, K1] B1;
  matrix[2, K1] B_test1;
  matrix[N, K2] B2;
  matrix[2, K2] B_test2;
  real<lower=0> var_inflate;
  real<lower=0> a0;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K1] b_spline;
  vector[2+K1] dbeta_raw;
  vector<lower=0>[2+K1] lambda_dbeta;
  real<lower=0> tau_dbeta;
  
  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  matrix[n_patients, K2] r_spline_raw;
  vector[n_patients] r0_rem_raw;
  vector[n_patients] r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem_raw;
  
  vector<lower=0>[n_patients] lambda_dr;
  real<lower=0> tau_dr;

  real<lower=0.8, upper=6> sigma_ref;
  real<lower=0.8, upper=6> sigma_rem;
  vector<lower=0, upper=80>[n_cls-2] sigma_d_raw;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  real<lower=0> sigma_r_spline;
  
  // mixture parameters
  real<lower=0> e0;
  simplex[n_cls] pi; 
}

transformed parameters {
  vector[n_cls-2] sigma_d = sort_asc(sigma_d_raw);
  vector<lower=0>[n_cls] sigma; //remote data
  vector[N] mu_y;
  vector[2+K1] delta_beta = dbeta_raw .* lambda_dbeta * tau_dbeta;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  matrix[n_patients, K2] r_spline = sigma_r_spline * r_spline_raw;
  
  vector[n_patients] sigma_r0_rem = sqrt(square(sigma_r0) + square(lambda_dr) * square(tau_dr));
  vector[n_patients] sigma_r1_rem = sqrt(square(sigma_r1) + square(lambda_dr) * square(tau_dr));
  matrix[n_patients, K2] sigma_r_spline_rem = sqrt(square(sigma_r_spline) + square(rep_matrix(lambda_dr, K2)) .* square(tau_dr));
  
  vector[n_patients] r0_rem = sigma_r0_rem .* r0_rem_raw;
  vector[n_patients] r1_rem = sigma_r1_rem .* r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem  = sigma_r_spline_rem .* r_spline_rem_raw;
  
  vector[N] mu_fixed = beta0 + beta_bl * x_bl + beta1 * time + B1 * b_spline;
  
  sigma[1] = sigma_rem;
  sigma[2] = sqrt(var_inflate) * sigma_ref;
  sigma[3:n_cls] = sqrt(square(sigma[2]) + square(sigma_d));

  
  for (i in idx_ref) {
    mu_y[i] = mu_fixed[i] + r0[patient_id[i]] + r1[patient_id[i]] * time[i] + dot_product(r_spline[patient_id[i], ], B2[i, ]);
  }
  for (i in idx_rem) {
    mu_y[i] = mu_fixed[i] + delta_beta[1] + delta_beta[2] * time[i] + B1[i, ] * delta_beta[3:(2+K1)] + r0_rem[patient_id[i]] + r1_rem[patient_id[i]] * time[i] + dot_product(r_spline_rem[patient_id[i], ], B2[i, ]);
  }
}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
  sigma_rem ~ student_t(1, 0, 5); //student_t(1, 0, 5);
  sigma_d_raw ~ student_t(3, 0, 5);
  
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,)
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  sigma_r_spline ~ student_t(1, 0, 5);
  
  dbeta_raw ~ std_normal();
  lambda_dbeta ~ lambda_prior;
  tau_dbeta ~ tau_prior;
  
  r0_rem_raw ~ std_normal();
  r1_rem_raw ~ std_normal();
  to_vector(r_spline_rem_raw) ~ std_normal();
  
  lambda_dr ~ lambda_prior;
  tau_dr ~ tau_prior;
  
  // mixture weight
  e0 ~ gamma(a0, a0 * n_cls);
  pi ~ dirichlet(rep_vector(e0, n_cls));

  // likelihood
  //reference
  y[idx_ref] ~ normal(mu_y[idx_ref], sigma_ref);
  
  // remote: GMM
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  array[N] int<lower=1, upper=n_cls> z;

  z[idx_ref] = rep_array(1, size(idx_ref));
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    z[i] = categorical_rng(softmax(lps));
  }
  
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  vector[K2] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (k in 1:K2) {
      r_spline_test[k] = normal_rng(0, sigma_r_spline);
    }
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test1[2, ], b_spline);
    real r_spline_contrib = dot_product(r_spline_test, B_test2[2, ]);
    
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12 + r_spline_contrib;
    y_pred_t12[i] = mu_y_test[i];
  }
  
  
}

"


### sGMM.transfer: v6: ri, ri*, ridge reg ============
sGMM.transfer6.semisep4 <- "
data {
  int<lower=0> N;
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K1;
  int<lower=0> K2;
  int<lower=2> n_cls;  
  int<lower=0> n_ref;
  int<lower=0> n_rem;
  vector[N] y;     
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;  
  array[n_ref] int idx_ref;  
  array[n_rem] int idx_rem;  
  matrix[N, K1] B1;
  matrix[2, K1] B_test1;
  matrix[N, K2] B2;
  matrix[2, K2] B_test2;
  real<lower=0> var_inflate;
  real<lower=0> a0;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K1] b_spline;
  vector[2+K1] dbeta_raw;
  vector<lower=0>[2+K1] lambda_dbeta;
  real<lower=0> tau_dbeta;
  
  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  matrix[n_patients, K2] r_spline_raw;
  vector[n_patients] r0_rem_raw;
  vector[n_patients] r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem_raw;

  real<lower=0.8, upper=6> sigma_ref;
  real<lower=0.8, upper=6> sigma_rem;
  vector<lower=0, upper=80>[n_cls-2] sigma_d_raw;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  real<lower=0> sigma_r_spline;
  real<lower=0> sigma_dr0;
  real<lower=0> sigma_dr1;
  real<lower=0> sigma_dr_spline;
  
  // mixture parameters
  real<lower=0> e0;
  simplex[n_cls] pi; 
}

transformed parameters {
  vector[n_cls-2] sigma_d = sort_asc(sigma_d_raw);
  vector<lower=0>[n_cls] sigma; //remote data
  vector[N] mu_y;
  vector[2+K1] delta_beta = dbeta_raw .* lambda_dbeta * tau_dbeta;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  matrix[n_patients, K2] r_spline = sigma_r_spline * r_spline_raw;
  
  real<lower=0> sigma_r0_rem = sqrt(square(sigma_r0) + square(sigma_dr0));
  real<lower=0> sigma_r1_rem = sqrt(square(sigma_r1) + square(sigma_dr1));
  real<lower=0> sigma_r_spline_rem = sqrt(square(sigma_r_spline) + square(sigma_dr_spline));
  
  vector[n_patients] r0_rem = sigma_r0_rem * r0_rem_raw;
  vector[n_patients] r1_rem = sigma_r1_rem * r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem  = sigma_r_spline_rem * r_spline_rem_raw;
  
  vector[N] mu_fixed = beta0 + beta_bl * x_bl + beta1 * time + B1 * b_spline;
  
  sigma[1] = sigma_rem;
  sigma[2] = sqrt(var_inflate) * sigma_ref;
  sigma[3:n_cls] = sqrt(square(sigma[2]) + square(sigma_d));

  
  for (i in idx_ref) {
    mu_y[i] = mu_fixed[i] + r0[patient_id[i]] + r1[patient_id[i]] * time[i] + dot_product(r_spline[patient_id[i], ], B2[i, ]);
  }
  for (i in idx_rem) {
    mu_y[i] = mu_fixed[i] + delta_beta[1] + delta_beta[2] * time[i] + B1[i, ] * delta_beta[3:(2+K1)] + r0_rem[patient_id[i]] + r1_rem[patient_id[i]] * time[i] + dot_product(r_spline_rem[patient_id[i], ], B2[i, ]);
  }
}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
  sigma_rem ~ student_t(1, 0, 5); //student_t(1, 0, 5);
  sigma_d_raw ~ student_t(3, 0, 5);
  
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,)
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  sigma_r_spline ~ student_t(1, 0, 5);
  
  dbeta_raw ~ std_normal();
  lambda_dbeta ~ lambda_prior;
  tau_dbeta ~ tau_prior;
  
  r0_rem_raw ~ std_normal();
  r1_rem_raw ~ std_normal();
  to_vector(r_spline_rem_raw) ~ std_normal();
  sigma_dr0 ~ student_t(1, 0, 5);
  sigma_dr1 ~ student_t(1, 0, 5);
  sigma_dr_spline ~ student_t(1, 0, 5);
  
  // mixture weight
  e0 ~ gamma(a0, a0 * n_cls);
  pi ~ dirichlet(rep_vector(e0, n_cls));

  // likelihood
  //reference
  y[idx_ref] ~ normal(mu_y[idx_ref], sigma_ref);
  
  // remote: GMM
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  array[N] int<lower=1, upper=n_cls> z;

  z[idx_ref] = rep_array(1, size(idx_ref));
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    z[i] = categorical_rng(softmax(lps));
  }
  
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  vector[K2] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (k in 1:K2) {
      r_spline_test[k] = normal_rng(0, sigma_r_spline);
    }
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test1[2, ], b_spline);
    real r_spline_contrib = dot_product(r_spline_test, B_test2[2, ]);
    
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12 + r_spline_contrib;
    y_pred_t12[i] = mu_y_test[i];
  }
  
  
}

"

### sGMM.transfer: v6, ri, ri*, spike-and-slab ============
sGMM.transfer6.semisep5 <- "
data {
  int<lower=0> N;
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K1;
  int<lower=0> K2;
  int<lower=2> n_cls;
  int<lower=0> n_ref;
  int<lower=0> n_rem;
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  array[N_obs] int idx_obs;
  array[N_mis] int idx_mis;
  vector[N_obs] y_obs;
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;
  array[n_ref] int idx_ref;
  array[n_rem] int idx_rem;
  matrix[N, K1] B1;
  matrix[2, K1] B_test1;
  matrix[N, K2] B2;
  matrix[2, K2] B_test2;
  real<lower=0> var_inflate;
  real<lower=0> a0;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K1] b_spline;
  vector[2+K1] dbeta_raw;
  vector<lower=0>[2+K1] lambda_dbeta;
  real<lower=0> tau_dbeta;
  
  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  matrix[n_patients, K2] r_spline_raw;
  vector[n_patients] r0_rem_raw;
  vector[n_patients] r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem_raw;

  real<lower=0.8, upper=6> sigma_ref;
  real<lower=0.8, upper=6> sigma_rem;
  vector<lower=0, upper=80>[n_cls-2] sigma_d_raw;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  real<lower=0> sigma_r_spline;
  real<lower=0> prec_r0_rem;
  real<lower=0> prec_r1_rem;
  real<lower=0> prec_r_spline_rem;
  real<lower=-1, upper=1> rho_r0;
  real<lower=-1, upper=1> rho_r1;
  real<lower=-1, upper=1> rho_r_spline;
  real<lower=0, upper=1> w_r0;
  real<lower=0, upper=1> w_r1;
  real<lower=0, upper=1> w_r_spline;

  // mixture parameters
  real<lower=0> e0;
  simplex[n_cls] pi;
  vector[N_mis] y_mis;
}

transformed parameters {
  vector[N] y;
  vector[n_cls-2] sigma_d = sort_asc(sigma_d_raw);
  vector<lower=0>[n_cls] sigma; //remote data
  vector[N] mu_y;
  vector[2+K1] delta_beta = dbeta_raw .* lambda_dbeta * tau_dbeta;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  matrix[n_patients, K2] r_spline = sigma_r_spline * r_spline_raw;
  real<lower=0> sigma_r0_rem = 1.0 / sqrt(prec_r0_rem);
  real<lower=0> sigma_r1_rem = 1.0 / sqrt(prec_r1_rem);
  real<lower=0> sigma_r_spline_rem = 1.0 / sqrt(prec_r_spline_rem);
  vector[n_patients] r0_rem = sigma_r0_rem * (rho_r0 * r0_raw + sqrt(1 - square(rho_r0)) * r0_rem_raw);
  vector[n_patients] r1_rem = sigma_r1_rem * (rho_r1 * r1_raw + sqrt(1 - square(rho_r1)) * r1_rem_raw);
  matrix[n_patients, K2] r_spline_rem  = sigma_r_spline_rem * (rho_r_spline * r_spline_raw + sqrt(1 - square(rho_r_spline)) * r_spline_rem_raw);

  vector[N] mu_fixed = beta0 + beta_bl * x_bl + beta1 * time + B1 * b_spline;

  y[idx_obs] = y_obs;
  y[idx_mis] = y_mis;
  sigma[1] = sigma_rem;
  sigma[2] = sqrt(var_inflate) * sigma_ref;
  sigma[3:n_cls] = sqrt(square(sigma[2]) + square(sigma_d));

  
  for (i in idx_ref) {
    mu_y[i] = mu_fixed[i] + r0[patient_id[i]] + r1[patient_id[i]] .* time[i] + dot_product(r_spline[patient_id[i], ], B2[i, ]);
  }
  for (i in idx_rem) {
    mu_y[i] = mu_fixed[i] + delta_beta[1] + delta_beta[2] * time[i] + B1[i, ] * delta_beta[3:(2+K1)] + r0_rem[patient_id[i]] + r1_rem[patient_id[i]] .* time[i] + dot_product(r_spline_rem[patient_id[i], ], B2[i, ]);
  }
}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
  sigma_rem ~ student_t(1, 0, 5); //student_t(1, 0, 5);
  sigma_d_raw ~ student_t(3, 0, 5);
  
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,)
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  sigma_r_spline ~ student_t(1, 0, 5);
  
  dbeta_raw ~ std_normal();
  lambda_dbeta ~ lambda_prior;
  tau_dbeta ~ tau_prior;
  
  r0_rem_raw ~ std_normal();
  r1_rem_raw ~ std_normal();
  to_vector(r_spline_rem_raw) ~ std_normal();
  rho_r0 ~ normal(0, 1);
  rho_r1 ~ normal(0, 1);
  rho_r_spline ~ normal(0, 1);

  w_r0 ~ beta(1, 1);
  target += log_mix(w_r0,
                  gamma_lpdf(prec_r0_rem | 100, 100 * square(sigma_r0)),
                  gamma_lpdf(prec_r0_rem | a_slab, b_slab)
                 );
                 
  w_r1 ~ beta(1, 1);
  target += log_mix(w_r1,
                  gamma_lpdf(prec_r1_rem | 100, 100 * square(sigma_r1)),
                  gamma_lpdf(prec_r1_rem | a_slab, b_slab)
                 );
                 
  w_r_spline ~ beta(1, 1);
  target += log_mix(w_r_spline,
                  gamma_lpdf(prec_r_spline_rem | 100, 100 * square(sigma_r_spline)),
                  gamma_lpdf(prec_r_spline_rem | a_slab, b_slab)
                 );
  
  // mixture weight
  e0 ~ gamma(a0, a0 * n_cls);
  pi ~ dirichlet(rep_vector(e0, n_cls));

  // likelihood
  //reference
  y[idx_ref] ~ normal(mu_y[idx_ref], sigma_ref);
  
  // remote: GMM
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  array[N] int<lower=1, upper=n_cls> z;

  z[idx_ref] = rep_array(1, size(idx_ref));
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    z[i] = categorical_rng(softmax(lps));
  }
  
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  vector[K2] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (k in 1:K2) {
      r_spline_test[k] = normal_rng(0, sigma_r_spline);
    }
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test1[2, ], b_spline);
    real r_spline_contrib = dot_product(r_spline_test, B_test2[2, ]);
    
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12 + r_spline_contrib;
    y_pred_t12[i] = mu_y_test[i];
  }
  
  
}

"

### sGMM.transfer: v6, ri, ri*, spike-and-slab, no correlation ============
sGMM.transfer6.semisep5.nocor <- "
data {
  int<lower=0> N;
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K1;
  int<lower=0> K2;
  int<lower=2> n_cls;
  int<lower=0> n_ref;
  int<lower=0> n_rem;
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  array[N_obs] int idx_obs;
  array[N_mis] int idx_mis;
  vector[N_obs] y_obs;
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;
  array[n_ref] int idx_ref;
  array[n_rem] int idx_rem;
  matrix[N, K1] B1;
  matrix[2, K1] B_test1;
  matrix[N, K2] B2;
  matrix[2, K2] B_test2;
  real<lower=0> var_inflate;
  real<lower=0> a0;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K1] b_spline;
  vector[2+K1] dbeta_raw;
  vector<lower=0>[2+K1] lambda_dbeta;
  real<lower=0> tau_dbeta;

  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  matrix[n_patients, K2] r_spline_raw;
  vector[n_patients] r0_rem_raw;
  vector[n_patients] r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem_raw;

  real<lower=0.8, upper=6> sigma_ref;
  real<lower=0.8, upper=6> sigma_rem;
  vector<lower=0, upper=80>[n_cls-2] sigma_d_raw;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  real<lower=0> sigma_r_spline;
  real<lower=0> prec_r0_rem;
  real<lower=0> prec_r1_rem;
  real<lower=0> prec_r_spline_rem;
  real<lower=0, upper=1> w_r0;
  real<lower=0, upper=1> w_r1;
  real<lower=0, upper=1> w_r_spline;

  // mixture parameters
  real<lower=0> e0;
  simplex[n_cls] pi;
  vector[N_mis] y_mis;
}

transformed parameters {
  vector[N] y;
  vector[n_cls-2] sigma_d = sort_asc(sigma_d_raw);
  vector<lower=0>[n_cls] sigma; //remote data
  vector[N] mu_y;
  vector[2+K1] delta_beta = dbeta_raw .* lambda_dbeta * tau_dbeta;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  matrix[n_patients, K2] r_spline = sigma_r_spline * r_spline_raw;
  real<lower=0> sigma_r0_rem = 1.0 / sqrt(prec_r0_rem);
  real<lower=0> sigma_r1_rem = 1.0 / sqrt(prec_r1_rem);
  real<lower=0> sigma_r_spline_rem = 1.0 / sqrt(prec_r_spline_rem);
  vector[n_patients] r0_rem = sigma_r0_rem * r0_rem_raw;
  vector[n_patients] r1_rem = sigma_r1_rem * r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem = sigma_r_spline_rem * r_spline_rem_raw;

  vector[N] mu_fixed = beta0 + beta_bl * x_bl + beta1 * time + B1 * b_spline;

  y[idx_obs] = y_obs;
  y[idx_mis] = y_mis;
  sigma[1] = sigma_rem;
  sigma[2] = sqrt(var_inflate) * sigma_ref;
  sigma[3:n_cls] = sqrt(square(sigma[2]) + square(sigma_d));

  for (i in idx_ref) {
    mu_y[i] = mu_fixed[i] + r0[patient_id[i]] + r1[patient_id[i]] .* time[i] + dot_product(r_spline[patient_id[i], ], B2[i, ]);
  }
  for (i in idx_rem) {
    mu_y[i] = mu_fixed[i] + delta_beta[1] + delta_beta[2] * time[i] + B1[i, ] * delta_beta[3:(2+K1)] + r0_rem[patient_id[i]] + r1_rem[patient_id[i]] .* time[i] + dot_product(r_spline_rem[patient_id[i], ], B2[i, ]);
  }
}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
  sigma_rem ~ student_t(1, 0, 5);
  sigma_d_raw ~ student_t(3, 0, 5);

  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5);

  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  sigma_r_spline ~ student_t(1, 0, 5);

  dbeta_raw ~ std_normal();
  lambda_dbeta ~ lambda_prior;
  tau_dbeta ~ tau_prior;

  r0_rem_raw ~ std_normal();
  r1_rem_raw ~ std_normal();
  to_vector(r_spline_rem_raw) ~ std_normal();

  w_r0 ~ beta(1, 1);
  target += log_mix(w_r0,
                  gamma_lpdf(prec_r0_rem | 100, 100 * square(sigma_r0)),
                  gamma_lpdf(prec_r0_rem | a_slab, b_slab)
                 );

  w_r1 ~ beta(1, 1);
  target += log_mix(w_r1,
                  gamma_lpdf(prec_r1_rem | 100, 100 * square(sigma_r1)),
                  gamma_lpdf(prec_r1_rem | a_slab, b_slab)
                 );

  w_r_spline ~ beta(1, 1);
  target += log_mix(w_r_spline,
                  gamma_lpdf(prec_r_spline_rem | 100, 100 * square(sigma_r_spline)),
                  gamma_lpdf(prec_r_spline_rem | a_slab, b_slab)
                 );

  // mixture weight
  e0 ~ gamma(a0, a0 * n_cls);
  pi ~ dirichlet(rep_vector(e0, n_cls));

  // likelihood
  y[idx_ref] ~ normal(mu_y[idx_ref], sigma_ref);

  // remote: GMM
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  array[N] int<lower=1, upper=n_cls> z;

  z[idx_ref] = rep_array(1, size(idx_ref));
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    z[i] = categorical_rng(softmax(lps));
  }

  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  vector[K2] r_spline_test;

  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (k in 1:K2) {
    r_spline_test[k] = normal_rng(0, sigma_r_spline);
  }
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test1[2, ], b_spline);
    real r_spline_contrib = dot_product(r_spline_test, B_test2[2, ]);

    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12 + r_spline_contrib;
    y_pred_t12[i] = mu_y_test[i];
  }
}

"


### over-parameterization versions with delta: discarded

### sGMM.transfer: v6, hetero RE ============
sGMM.transfer6.heteroRE <- "
data {
  int<lower=0> N;
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K1;
  int<lower=0> K2;
  int<lower=2> n_cls;  
  int<lower=0> n_ref;
  int<lower=0> n_rem;
  vector[N] y;     
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;  
  array[n_ref] int idx_ref;  
  array[n_rem] int idx_rem;  
  matrix[N, K1] B1;
  matrix[2, K1] B_test1;
  matrix[N, K2] B2;
  matrix[2, K2] B_test2;
  real<lower=0> var_inflate;
  real<lower=0> a0;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K1] b_spline;
  vector[2+K1] dbeta_raw;
  vector<lower=0>[2+K1] lambda_dbeta;
  real<lower=0> tau_dbeta;
  
  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  matrix[n_patients, K2] r_spline_raw;
  vector[n_patients] r0_rem_raw;
  vector[n_patients] r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem_raw;
  vector[n_patients] dr0_raw;
  vector[n_patients] dr1_raw;
  matrix[n_patients, K2] delta_r_spline_raw;
  vector<lower=0>[n_patients] lambda_dr0;
  real<lower=0> tau_dr0;
  vector<lower=0>[n_patients] lambda_dr1;
  real<lower=0> tau_dr1;
  matrix<lower=0>[n_patients, K2] lambda_dr_spline;
  real<lower=0> tau_dr_spline;

  real<lower=0.8, upper=6> sigma_ref;
  real<lower=0.8, upper=6> sigma_rem;
  //positive_ordered[n_cls-2] sigma_d;
  vector<lower=0, upper=80>[n_cls-2] sigma_d_raw;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  real<lower=0> sigma_r_spline;
  real<lower=0> sigma_r0_rem;
  real<lower=0> sigma_r1_rem;
  real<lower=0> sigma_r_spline_rem;
  
  // mixture parameters
  real<lower=0> e0;
  simplex[n_cls] pi; 
}

transformed parameters {
  vector[n_cls-2] sigma_d = sort_asc(sigma_d_raw);
  vector<lower=0>[n_cls] sigma; //remote data
  vector[N] mu_y;
  vector[2+K1] delta_beta = dbeta_raw .* lambda_dbeta * tau_dbeta;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  vector[n_patients] r0_rem = sigma_r0_rem * r0_rem_raw;
  vector[n_patients] r1_rem = sigma_r1_rem * r1_rem_raw;
  matrix[n_patients, K2] r_spline = sigma_r_spline * r_spline_raw;
  matrix[n_patients, K2] r_spline_rem  = sigma_r_spline_rem * r_spline_rem_raw;
  
  vector[n_patients] delta_r0_rem = dr0_raw .* lambda_dr0 * tau_dr0;
  vector[n_patients] delta_r1_rem = dr1_raw .* lambda_dr1 * tau_dr1;
  matrix[n_patients, K2] delta_r_spline_rem = delta_r_spline_raw .* lambda_dr_spline * tau_dr_spline;
  vector[N] mu_fixed = beta0 + beta_bl * x_bl + beta1 * time + B1 * b_spline;
  
  sigma[1] = sigma_rem;
  sigma[2] = sqrt(var_inflate) * sigma_ref;
  sigma[3:n_cls] = sqrt(square(sigma[2]) + square(sigma_d));

  for (i in idx_ref) {
    mu_y[i] = mu_fixed[i] + r0[patient_id[i]] + r1[patient_id[i]] * time[i] + dot_product(r_spline[patient_id[i], ], B2[i, ]);
  }
  for (i in idx_rem) {
    mu_y[i] = mu_fixed[i] + delta_beta[1] + delta_beta[2] * time[i] + B1[i, ] * delta_beta[3:(2+K1)] + r0_rem[patient_id[i]] + delta_r0_rem[patient_id[i]] + r1_rem[patient_id[i]] * time[i] + delta_r1_rem[patient_id[i]] * time[i] + dot_product(r_spline_rem[patient_id[i], ], B2[i, ]) + dot_product(delta_r_spline_rem[patient_id[i], ], B2[i, ]);
  }

}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
  sigma_rem ~ sigma_prior; //student_t(1, 0, 5);
  sigma_d_raw ~ student_t(3, 0, 5);
  
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,)
  
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  sigma_r_spline ~ student_t(1, 0, 5);
  
  r0_rem_raw ~ std_normal();
  r1_rem_raw ~ std_normal();
  to_vector(r_spline_rem_raw) ~ std_normal();
  sigma_r0_rem ~ student_t(1, 0, 5);
  sigma_r1_rem ~ student_t(1, 0, 5);
  sigma_r_spline_rem ~ student_t(1, 0, 5);
  
  dbeta_raw ~ std_normal();
  lambda_dbeta ~ lambda_prior;
  tau_dbeta ~ tau_prior;
  
  dr0_raw ~ std_normal();
  lambda_dr0 ~ lambda_prior;
  tau_dr0 ~ tau_prior2;
  
  dr1_raw ~ std_normal();
  lambda_dr1 ~ lambda_prior;
  tau_dr1 ~ tau_prior2;
  
  to_vector(delta_r_spline_raw) ~ std_normal();
  to_vector(lambda_dr_spline) ~ lambda_prior;
  tau_dr_spline ~ tau_prior2;
  
  // mixture weight
  e0 ~ gamma(a0, a0 * n_cls);
  pi ~ dirichlet(rep_vector(e0, n_cls));

  // likelihood
  //reference
  y[idx_ref] ~ normal(mu_y[idx_ref], sigma_ref);
  
  // remote: GMM
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  array[N] int<lower=1, upper=n_cls> z;

  z[idx_ref] = rep_array(1, size(idx_ref));
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    z[i] = categorical_rng(softmax(lps));
  }
  
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  vector[K2] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (k in 1:K2) {
      r_spline_test[k] = normal_rng(0, sigma_r_spline);
    }
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test1[2, ], b_spline);
    real r_spline_contrib = dot_product(r_spline_test, B_test2[2, ]);
    
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12 + r_spline_contrib;
    y_pred_t12[i] = mu_y_test[i];
  }
  
  
}

"
### sGMM.transfer: v6, hetero RE, subject-level shrinkage ============
sGMM.transfer6.heteroRE2 <- "
data {
  int<lower=0> N;
  int<lower=0> N_test;
  int<lower=0> n_patients;
  int<lower=0> K1;
  int<lower=0> K2;
  int<lower=2> n_cls;  
  int<lower=0> n_ref;
  int<lower=0> n_rem;
  vector[N] y;     
  vector[N] x_bl;
  vector[N_test] x_bl_test;
  vector[N] time;
  array[N] int<lower=1, upper=n_patients> patient_id;  
  array[n_ref] int idx_ref;  
  array[n_rem] int idx_rem;  
  matrix[N, K1] B1;
  matrix[2, K1] B_test1;
  matrix[N, K2] B2;
  matrix[2, K2] B_test2;
  real<lower=0> var_inflate;
  real<lower=0> a0;
}


parameters {
  // fixed effect
  real beta0;
  real beta1;
  real beta_bl;
  vector[K1] b_spline;
  vector[2+K1] dbeta_raw;
  vector<lower=0>[2+K1] lambda_dbeta;
  real<lower=0> tau_dbeta;
  
  // random effect
  vector[n_patients] r0_raw;
  vector[n_patients] r1_raw;
  matrix[n_patients, K2] r_spline_raw;
  vector[n_patients] r0_rem_raw;
  vector[n_patients] r1_rem_raw;
  matrix[n_patients, K2] r_spline_rem_raw;
  vector[n_patients] dr0_raw;
  vector[n_patients] dr1_raw;
  matrix[n_patients, K2] delta_r_spline_raw;
  vector<lower=0>[n_patients] lambda_dr;
  real<lower=0> tau_dr;

  real<lower=0.8, upper=6> sigma_ref;
  real<lower=0.8, upper=6> sigma_rem;
  //positive_ordered[n_cls-2] sigma_d;
  vector<lower=0, upper=80>[n_cls-2] sigma_d_raw;
  real<lower=0> sigma_b_spline;
  real<lower=0> sigma_r0;
  real<lower=0> sigma_r1;
  real<lower=0> sigma_r_spline;
  real<lower=0> sigma_r0_rem;
  real<lower=0> sigma_r1_rem;
  real<lower=0> sigma_r_spline_rem;
  
  // mixture parameters
  real<lower=0> e0;
  simplex[n_cls] pi; 
}

transformed parameters {
  vector[n_cls-2] sigma_d = sort_asc(sigma_d_raw);
  vector<lower=0>[n_cls] sigma; //remote data
  vector[N] mu_y;
  vector[2+K1] delta_beta = dbeta_raw .* lambda_dbeta * tau_dbeta;
  vector[n_patients] r0 = sigma_r0 * r0_raw;
  vector[n_patients] r1 = sigma_r1 * r1_raw;
  vector[n_patients] r0_rem = sigma_r0_rem * r0_rem_raw;
  vector[n_patients] r1_rem = sigma_r1_rem * r1_rem_raw;
  matrix[n_patients, K2] r_spline = sigma_r_spline * r_spline_raw;
  matrix[n_patients, K2] r_spline_rem  = sigma_r_spline_rem * r_spline_rem_raw;
  
  vector[n_patients] delta_r0_rem = dr0_raw .* lambda_dr * tau_dr;
  vector[n_patients] delta_r1_rem = dr1_raw .* lambda_dr * tau_dr;
  matrix[n_patients, K2] delta_r_spline_rem = delta_r_spline_raw .* rep_matrix(lambda_dr, K2) * tau_dr;
  vector[N] mu_fixed = beta0 + beta_bl * x_bl + beta1 * time + B1 * b_spline;
  
  sigma[1] = sigma_rem;
  sigma[2] = sqrt(var_inflate) * sigma_ref;
  sigma[3:n_cls] = sqrt(square(sigma[2]) + square(sigma_d));

  for (i in idx_ref) {
    mu_y[i] = mu_fixed[i] + r0[patient_id[i]] + r1[patient_id[i]] * time[i] + dot_product(r_spline[patient_id[i], ], B2[i, ]);
  }
  for (i in idx_rem) {
    mu_y[i] = mu_fixed[i] + delta_beta[1] + delta_beta[2] * time[i] + B1[i, ] * delta_beta[3:(2+K1)] + r0_rem[patient_id[i]] + delta_r0_rem[patient_id[i]] + r1_rem[patient_id[i]] * time[i] + delta_r1_rem[patient_id[i]] * time[i] + dot_product(r_spline_rem[patient_id[i], ], B2[i, ]) + dot_product(delta_r_spline_rem[patient_id[i], ], B2[i, ]);
  }

}

model {
  // mixture parameters
  sigma_ref ~ student_t(1, 0, 5);
  sigma_rem ~ sigma_prior; //student_t(1, 0, 5);
  sigma_d_raw ~ student_t(3, 0, 5);
  
  //fixed coef
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta_bl ~ normal(0, 10);
  b_spline ~ normal(0, sigma_b_spline);
  sigma_b_spline ~ student_t(1, 0, 5); //jags: dt(0, 5^(-2), 1)T(0,)
  
  
  //random coef
  r0_raw ~ std_normal();
  r1_raw ~ std_normal();
  to_vector(r_spline_raw) ~ std_normal();
  sigma_r0 ~ student_t(1, 0, 5);
  sigma_r1 ~ student_t(1, 0, 5);
  sigma_r_spline ~ student_t(1, 0, 5);
  
  r0_rem_raw ~ std_normal();
  r1_rem_raw ~ std_normal();
  to_vector(r_spline_rem_raw) ~ std_normal();
  sigma_r0_rem ~ student_t(1, 0, 5);
  sigma_r1_rem ~ student_t(1, 0, 5);
  sigma_r_spline_rem ~ student_t(1, 0, 5);
  
  dbeta_raw ~ std_normal();
  lambda_dbeta ~ lambda_prior;
  tau_dbeta ~ tau_prior;
  
  dr0_raw ~ std_normal();
  dr1_raw ~ std_normal();
  to_vector(delta_r_spline_raw) ~ std_normal();
  lambda_dr ~ lambda_prior;
  tau_dr ~ tau_prior2;
  
  // mixture weight
  e0 ~ gamma(a0, a0 * n_cls);
  pi ~ dirichlet(rep_vector(e0, n_cls));

  // likelihood
  //reference
  y[idx_ref] ~ normal(mu_y[idx_ref], sigma_ref);
  
  // remote: GMM
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  array[N] int<lower=1, upper=n_cls> z;

  z[idx_ref] = rep_array(1, size(idx_ref));
  for (i in idx_rem) {
    vector[n_cls] lps = log(pi);
    for (cc in 1:n_cls) {
      lps[cc] += normal_lpdf(y[i] | mu_y[i], sigma[cc]);
    }
    z[i] = categorical_rng(softmax(lps));
  }
  
  array[N_test] real y_pred_t12;
  array[N_test] real mu_fixed_test;
  array[N_test] real mu_y_test;
  real r0_test;
  real r1_test;
  vector[K2] r_spline_test;
  
  r0_test = normal_rng(0, sigma_r0);
  r1_test = normal_rng(0, sigma_r1);
  for (k in 1:K2) {
      r_spline_test[k] = normal_rng(0, sigma_r_spline);
    }
  for (i in 1:N_test) {
    mu_fixed_test[i] = beta0 + beta_bl * x_bl_test[i] + beta1 * 12 + dot_product(B_test1[2, ], b_spline);
    real r_spline_contrib = dot_product(r_spline_test, B_test2[2, ]);
    
    mu_y_test[i] = mu_fixed_test[i] + r0_test + r1_test * 12 + r_spline_contrib;
    y_pred_t12[i] = mu_y_test[i];
  }
  
  
}

"




