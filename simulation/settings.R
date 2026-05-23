beta_delta_setting <- matrix(
  c(0, 0, 0, 
    -2, 0, 0, 
    0, -2, 0), byrow = TRUE, nrow = 3
)
all.config <- bind_rows(
  # consistent scenario under type 3
  expand.grid(
    model_type = c(1:5),
    gmm_type = 1,
    n_patients = 80,
    beta_delta_idx = 1, 
    rho_source = c(0.5, 0.8), 
    re_factor = c(1), 
    remote_inflation = 0,
    outlier_rate1 = 0.2, # deprecated
    outlier_rate2 = 0.5, # deprecated
    outlier_type = 3,
    error_type = 1, # iid
    uni_modal = TRUE, # deprecated
    high_freq = FALSE,
    DGM_type = 1
  ),
  # gmm
  expand.grid(
    model_type = 6,
    gmm_type = c(3:5),
    n_patients = 80,
    beta_delta_idx = 1, 
    rho_source = c(0.5, 0.8), 
    re_factor = c(1), 
    remote_inflation = 0,
    outlier_rate1 = 0.2,
    outlier_rate2 = 0.5,
    outlier_type = 3,
    error_type = 1, # iid
    uni_modal = TRUE,
    high_freq = FALSE,
    DGM_type = 1
  ),
  # outlier 1: sporadic
  # rlmm
  expand.grid(
    model_type = c(1:5),
    gmm_type = 1,
    n_patients = 80,
    beta_delta_idx = 1:3, 
    rho_source = c(0.5, 0.8), 
    re_factor = c(1, 2), 
    remote_inflation = 0,
    outlier_rate1 = 0.5,
    outlier_rate2 = 0.2,
    outlier_type = 1,
    error_type = 1, # iid
    uni_modal = TRUE,
    high_freq = FALSE,
    DGM_type = 1
  ),
  # gmm
  expand.grid(
    model_type = 6,
    gmm_type = c(3:5),
    n_patients = 80,
    beta_delta_idx = 1:3, 
    rho_source = c(0.5, 0.8), 
    re_factor = c(1, 2), 
    remote_inflation = 0,
    outlier_rate1 = 0.5,
    outlier_rate2 = 0.2,
    outlier_type = 1,
    error_type = 1, # iid
    uni_modal = TRUE,
    high_freq = FALSE,
    DGM_type = 1
  ),
  # outlier 2: uni-modal recurring
  # rlmm
  expand.grid(
    model_type = c(1:5),
    gmm_type = 1,
    n_patients = 80,
    beta_delta_idx = 1:3, 
    rho_source = c(0.5, 0.8), 
    re_factor = c(1, 2), 
    remote_inflation = 0,
    outlier_rate1 = 0.2,
    outlier_rate2 = 0.5,
    outlier_type = 2,
    error_type = 1, # iid
    uni_modal = c(TRUE, FALSE),
    high_freq = c(FALSE),# c(TRUE, FALSE),
    DGM_type = 1
  ),
  # gmm
  expand.grid(
    model_type = 6,
    gmm_type = c(3:5),
    n_patients = 80,
    beta_delta_idx = 1:3, 
    rho_source = c(0.5, 0.8), 
    re_factor = c(1, 2), 
    remote_inflation = 0,
    outlier_rate1 = 0.2,
    outlier_rate2 = 0.5,
    outlier_type = 2,
    error_type = 1, # iid
    uni_modal = c(TRUE, FALSE),
    high_freq = c(FALSE),# c(TRUE, FALSE),
    DGM_type = 1
  ),
  # outlier 3: mixed outliers
  # rlmm
  expand.grid(
    model_type = c(1:5),
    gmm_type = 1,
    n_patients = 80,
    beta_delta_idx = 1:3, 
    rho_source = c(0.5, 0.8), 
    re_factor = c(1, 2), 
    remote_inflation = 0,
    outlier_rate1 = 0.2, # deprecated
    outlier_rate2 = 0.5, # deprecated
    outlier_type = 3,
    error_type = 1, # iid
    uni_modal = TRUE, # deprecated
    high_freq = c(FALSE),# c(TRUE, FALSE),
    DGM_type = 1
  ),
  # gmm
  expand.grid(
    model_type = 6,
    gmm_type = c(3:5),
    n_patients = 80,
    beta_delta_idx = 1:3, 
    rho_source = c(0.5, 0.8), 
    re_factor = c(1, 2), 
    remote_inflation = 0,
    outlier_rate1 = 0.2,
    outlier_rate2 = 0.5,
    outlier_type = 3,
    error_type = 1, # iid
    uni_modal = TRUE,
    high_freq = c(FALSE),# c(TRUE, FALSE),
    DGM_type = 1
  )
)

# which(all.config$gmm_type %in% c(3, 4, 5) & all.config$high_freq)
# which(all.config$gmm_type == 1)
