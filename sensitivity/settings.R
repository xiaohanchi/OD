beta_delta_setting <- matrix(
  c(0, 0, 0, 
    -2, 0, 0, 
    0, -2, 0), byrow = TRUE, nrow = 3
)
all.config <- bind_rows(
  # outlier 1: sporadic
  # gmm
  expand.grid(
    model_type = 6,
    gmm_type = 5,
    n_patients = 80,
    beta_delta_idx = c(1, 3), 
    rho_source = 0.8, # c(0.5, 0.8), 
    re_factor = 2, # c(1, 2), 
    remote_inflation = 0,
    outlier_rate1 = 0.5,
    outlier_rate2 = 0.2,
    outlier_type = 1,
    error_type = 1, # iid
    uni_modal = TRUE,
    high_freq = FALSE,
    gmm_ncls = c(5, 7, 10), 
    DGM_type = 1
  ),
  # outlier 2: uni-modal recurring
  # gmm
  expand.grid(
    model_type = 6,
    gmm_type = 5,
    n_patients = 80,
    beta_delta_idx = c(1, 3), 
    rho_source = 0.8, # c(0.5, 0.8), 
    re_factor = 2, # c(1, 2), 
    remote_inflation = 0,
    outlier_rate1 = 0.2,
    outlier_rate2 = 0.5,
    outlier_type = 2,
    error_type = 1, # iid
    uni_modal = c(TRUE, FALSE),
    high_freq = c(TRUE, FALSE),
    gmm_ncls = c(5, 7, 10), 
    DGM_type = 1
  ),
  # outlier 3: mixed outliers
  # gmm
  expand.grid(
    model_type = 6,
    gmm_type = 5,
    n_patients = 80,
    beta_delta_idx = c(1, 3), 
    rho_source = 0.8, # c(0.5, 0.8), 
    re_factor = 2, # c(1, 2), 
    remote_inflation = 0,
    outlier_rate1 = 0.2,
    outlier_rate2 = 0.5,
    outlier_type = 3,
    error_type = 1, # iid
    uni_modal = TRUE,
    high_freq = c(TRUE, FALSE),
    gmm_ncls = c(5, 7, 10), 
    DGM_type = 1
  )
)

# which(all.config$gmm_type %in% c(3, 4, 5) & all.config$high_freq)
# which(all.config$gmm_type == 1)
# which(all.config$model_type == 4)