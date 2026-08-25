rm(list = ls())

library(pacman)
p_load(
  coda, dplyr, emmeans, filelock, nlme, lme4, MASS, readr, mfp,
  robustlmm, rjags, runjags, rstan, splines, stringr
)
select <- dplyr::select

source("../stan_functions.R")
source("../settings.R")
source("../functions.R")

main_func <- function(group, n_patients, n_measure_o, n_measure_r, beta_delta, 
                      rho_source, re_factor, shift_rate, remote_inflation,
                      outlier_rate1, outlier_rate2, outlier_type,
                      uni_modal = TRUE, high_freq = FALSE,
                      DGM_type, model_type, gmm_ncls, gmm_type, error_type, 
                      seed = 123, rep = 1) {
  # model_type: 1 for robust smm; 2 for gmm and variants
  if (group == "trt") {
    sd_r <- c(2.5, 1)
    if (DGM_type == 1) beta <- c(0.20, -4.0, 0.10) else if (DGM_type == 2) beta <- c(-16.5, -30, -17)
  } else if (group == "ctrl") {
    beta <- c(0.10, -0.7, 0.05)
    sd_r <- c(1.8, 0.6)
  }
  
  if (error_type == 1) {
    corr_r <- 0.3
    rho_ar1 <- NA
    # sd_eps = c(4, sqrt(20))
    sd_eps = c(4, 5)
  } else if (error_type == 2) {
    corr_r <- 0.1
    # sd_r <- c(1, 1) # replaced in AR error
    rho_ar1 <- c(0.50, 0.40)
    sd_eps = c(4, 5)
    # sd_eps = c(3, sqrt(13))
  }

  sim_data <- generate_data(
    n_patients = n_patients, 
    n_measure = c("Onsite" = n_measure_o, "Remote" = n_measure_r),
    beta = beta, beta_delta = beta_delta, 
    sd_r = sd_r, corr_r = corr_r, rho_source = rho_source, re_factor = re_factor, 
    rho_ar1 = rho_ar1, sd_eps = sd_eps, shift_rate = shift_rate, remote_inflation = remote_inflation, 
    outlier_rate1 = outlier_rate1, outlier_rate2 = outlier_rate2,
    equal_remote_N = FALSE, outlier_type = outlier_type, 
    uni_modal = uni_modal, high_freq = high_freq, DGM_type = DGM_type, error_type = error_type, 
    seed = seed
  )
  
  oracle_data <- generate_data(
    n_patients = n_patients, 
    n_measure = c("Onsite" = n_measure_o, "Remote" = n_measure_r),
    beta = beta, beta_delta = beta_delta, 
    sd_r = sd_r, corr_r = corr_r, rho_source = rho_source, re_factor = re_factor, 
    rho_ar1 = rho_ar1, sd_eps = rep(sd_eps[1], 2), shift_rate = shift_rate, remote_inflation = 0, 
    outlier_rate1 = outlier_rate1, outlier_rate2 = outlier_rate2,
    equal_remote_N = FALSE, outlier_type = outlier_type, 
    uni_modal = uni_modal, high_freq = high_freq, DGM_type = DGM_type, error_type = error_type, 
    oracle = TRUE, seed = seed
  )
  test_data <- generate_data(
    n_patients = n_patients, 
    n_measure = c("Onsite" = n_measure_o, "Remote" = n_measure_r),
    beta = beta, beta_delta = beta_delta, 
    sd_r = sd_r, corr_r = corr_r, rho_source = rho_source, re_factor = re_factor, 
    rho_ar1 = rho_ar1, sd_eps = sd_eps, shift_rate = shift_rate, remote_inflation = remote_inflation, 
    outlier_rate1 = outlier_rate1, outlier_rate2 = outlier_rate2,
    equal_remote_N = FALSE, outlier_type = outlier_type, 
    uni_modal = uni_modal, high_freq = high_freq, DGM_type = DGM_type, error_type = error_type, 
    seed = 42
  ) %>% filter(source == "Onsite")


  if (model_type %in% c(1:5)) {
    convergence <- rep(-1, 3)
    if (model_type == 1) {
      # Robust SMM
      method <- "rsmm"
      res <- run_robustlmm(data_all = sim_data)
      y_fitted_all <- fitted(res$model_fit)
    } else if (model_type == 2) {
      method <- "onsite_smm"
      res <- run_lmm(data_all = sim_data, type = 1)
      y_fitted_all <- numeric(nrow(sim_data))
      # to match vector length
      y_fitted_all[sim_data$source == "Onsite"] <- fitted(res$model_fit)
    } else if (model_type == 3) {
      method <- "pool_smm"
      res <- run_lmm(data_all = sim_data, type = 2)
      y_fitted_all <- fitted(res$model_fit)
    } else if (model_type == 4) {
      method <- "oracle_smm"
      res <- run_lmm(data_all = oracle_data, type = 3, error_type = error_type)
      y_fitted_all <- fitted(res$model_fit)
    } else if (model_type == 5) {
      method <- "2stage_smm"
      res <- twostage_lmm(data_all = sim_data)
      y_fitted_all <- numeric(nrow(sim_data))
      y_fitted_all[!res$data_fit$outlier_est] <- fitted(res$model_fit)
    }
    y_pred <- predict(res$model_fit, newdata = test_data)
    
    if (model_type == 5) {
      muhat_fixed <- bind_cols(
        source =  res$data_fit$source[!res$data_fit$outlier_est], 
        time = res$data_fit$time[!res$data_fit$outlier_est], 
        mu_fixed = predict(res$model_fit, re.form = NA) %>% round(10)
      ) %>% filter(source == "Onsite") %>% select(-c(source)) %>% 
        distinct() # get (population level) fixed effect 
    } else if (model_type == 4) {
      if(error_type == 1) {
        muhat_fixed <- bind_cols(
          source =  res$data_fit$source, 
          time = res$data_fit$time, 
          mu_fixed = predict(res$model_fit, re.form = NA) %>% round(10)
        ) %>% filter(time %in% c(1:n_measure_o)) %>% 
          select(-c(source)) %>% distinct() 
      } else {
        muhat_fixed <- bind_cols(
          source =  res$data_fit$source, 
          time = res$data_fit$time, 
          mu_fixed = predict(res$model_fit, level = 0) %>% round(10)
        ) %>% filter(time %in% c(1:n_measure_o)) %>% 
          select(-c(source)) %>% distinct() 
      }
      
    } else {
      muhat_fixed <- bind_cols(
        source =  res$data_fit$source, 
        time = res$data_fit$time, 
        mu_fixed = predict(res$model_fit, re.form = NA) %>% round(10)
      ) %>% filter(source == "Onsite") %>% select(-c(source)) %>% 
        distinct() # get (population level) fixed effect 
    }
    ### CI 
    CI_res <- emmeans(res$model_fit, ~ time, 
            at = list(time = c(1:n_measure_o)),
            params = "k", data = res$data_fit, 
            level = 0.95, df = Inf) %>% confint()
    CI_length <- mean(CI_res$asymp.UCL - CI_res$asymp.LCL)
    CI_true <- (unique(sim_data[sim_data$source == "Onsite", c("time", "y_true_fixed")]) %>% arrange(time))$y_true_fixed
    CI_cov <- mean(CI_res$asymp.LCL <= CI_true & CI_res$asymp.UCL >= CI_true)
    
    CI_res_fwc <- emmeans(res$model_fit, ~ time, 
                          at = list(time = c(n_measure_o)),
                          params = "k", data = res$data_fit, level = 0.95, df = Inf) %>% 
      confint()
    CI_true_fwc <- unique(sim_data[sim_data$source == "Onsite", c("time", "y_true_fixed")]) %>% 
      arrange(time) %>% 
      filter(time == n_measure_o) %>%
      pull(y_true_fixed)
    CI_length_fwc <- CI_res_fwc$asymp.UCL - CI_res_fwc$asymp.LCL
    CI_cov_fwc <- (CI_res_fwc$asymp.LCL <= CI_true_fwc & CI_res_fwc$asymp.UCL >= CI_true_fwc) * 1
    HPD_length <- HPD_cov <- HPD_length_fwc <- HPD_cov_fwc <- NA
    
  } else if (model_type ==  6) {
    # GMM
    method <- paste0("gmm", gmm_type)
    res <- run_gmm(
      data_all = sim_data, test_data = test_data, n_cls = gmm_ncls,
      var_inflate = 5, a0 = 10, k1 = 5, type = gmm_type
    )
    muhat_fixed <- bind_cols(
      time = unique(res$data_fit$time[res$data_fit$group == "Onsite"]), 
      mu_fixed = colMeans(res$mu_fixed)
    ) 
    y_fitted_all <- res$mu_all
    y_pred <- res$y_pred
    # convergence rate
    convergence <- summary(res$stan_fit)$summary[c("beta0", "beta1", "sigma[1]"), c("n_eff", "Rhat")] %>% 
      round(., 2) %>% apply(2, max)
    convergence["converge"] <- ifelse(convergence["Rhat"] > 1.03, 0, 1)
    res$data_fit <- res$data_fit %>% rename(patient = id, source = group, y = y_obs)
    
    ### CI
    CI_res <- apply(res$mu_fixed, 2, quantile, probs = c(0.025, 0.975))
    CI_true <- (unique(sim_data[sim_data$source == "Onsite", c("time", "y_true_fixed")]) %>% arrange(time))$y_true_fixed
    CI_length <- mean(CI_res[2, ] - CI_res[1, ])
    CI_cov <- mean(CI_res[1, ] <= CI_true & CI_res[2, ] >= CI_true)
    HPD_res <- apply(res$mu_fixed, 2, function(x) {
      HPDinterval(as.mcmc(x), prob = 0.95)
    })
    HPD_length <- mean(HPD_res[2, ] - HPD_res[1, ])
    HPD_cov <- mean(HPD_res[1, ] <= CI_true & HPD_res[2, ] >= CI_true)
    
    CI_res_fwc <- (res$mu_fixed[, n_measure_o]) %>% quantile(probs = c(0.025, 0.975))
    CI_true_fwc <- unique(sim_data[sim_data$source == "Onsite", c("time", "y_true_fixed")]) %>% 
      arrange(time) %>% 
      filter(time == n_measure_o) %>%
      pull(y_true_fixed)
    CI_length_fwc <- diff(CI_res_fwc)
    CI_cov_fwc <- (CI_res_fwc[1] <= CI_true_fwc & CI_res_fwc[2] >= CI_true_fwc) * 1
    HPD_res_fwc <- HPDinterval(as.mcmc(res$mu_fixed[, n_measure_o]), prob = 0.95) 
    HPD_length_fwc <- HPD_res_fwc[1, "upper"] - HPD_res_fwc[1, "lower"]
    HPD_cov_fwc <- (HPD_res_fwc[1, "lower"] <= CI_true_fwc & HPD_res_fwc[1, "upper"] >= CI_true_fwc) * 1
    
  }
  
  sim_data <- sim_data %>% mutate(y_fitted_all = y_fitted_all) # to save
  if (model_type == 4) oracle_data <- oracle_data %>% mutate(y_fitted_all = y_fitted_all) # to save
  test_data <- test_data %>% mutate(y_pred = y_pred) 
  
  patient_data <-  res$data_fit %>% group_by(patient) %>% 
    summarise(
      eps_flag_3sigma = any(eps_flag_3sigma & outlier_est & source == "Remote"),
      eps_flag_2sigma = any(eps_flag_2sigma & outlier_est & source == "Remote"),
      outlier = any(outlier1 | outlier2), 
      outlier_det_3sigma = any((outlier1 | outlier2) & (!outlier_undetect_3sigma)), 
      outlier_det_2sigma = any((outlier1 | outlier2) & (!outlier_undetect_2sigma)), 
      outlier_est = any(outlier_est)
    )
  
  # point-wise outlier detection acc (FPR & FNR)
  fpr <- res$data_fit %>%
    summarise(fpr = mean(outlier_est[!outlier1 & !outlier2])) %>%
    pull()
  fnr <- res$data_fit %>%
    summarise(fnr = 1 - mean(outlier_est[outlier1 | outlier2])) %>%
    pull()
  cfpr_3sigma <- res$data_fit %>%
    summarise(cfpr = sum(outlier_est & !outlier1 & !outlier2 & !eps_flag_3sigma) / 
                sum(!outlier1 & !outlier2)) %>% pull()
  cfpr_2sigma <- res$data_fit %>%
    summarise(cfpr = sum(outlier_est & !outlier1 & !outlier2 & !eps_flag_2sigma) / 
                sum(!outlier1 & !outlier2)) %>% pull()
  cfnr_3sigma <- res$data_fit %>%
    summarise(cfnr = sum(!outlier_est & (outlier1 | outlier2) & !outlier_undetect_3sigma) / 
                sum(outlier1 | outlier2)) %>% pull()
  cfnr_2sigma <- res$data_fit %>%
    summarise(cfnr = sum(!outlier_est & (outlier1 | outlier2) & !outlier_undetect_2sigma) / 
                sum(outlier1 | outlier2)) %>% pull()
  
  # patient-wise outlier detection acc (FPR & FNR)
  fpr_p <- with(patient_data, {
    subset <- outlier_est[!outlier]
    c(sum(subset), length(subset))
  })
  fnr_p <- with(patient_data, {
    subset <- outlier_est[outlier]
    c(sum(!subset), length(subset))
  })
  error_p <- fpr_p + fnr_p
  
  cfpr_p <- with(patient_data, {
    subset <- outlier_est[!outlier]
    csubset <- outlier_est[(!outlier)&(!eps_flag_3sigma)]
    c(sum(csubset), length(subset))
  })
  cfnr_p <- with(patient_data, {
    subset <- outlier_est[outlier]
    csubset <- outlier_est[outlier_det_3sigma]
    c(sum(!csubset), length(subset))
  })
  cerror_3sigma_p <- cfpr_p + cfnr_p
  
  cerror_2sigma_p <- with(patient_data, {
    subset <- outlier_est[!outlier]
    csubset <- outlier_est[(!outlier)&(!eps_flag_2sigma)]
    c(sum(csubset), length(subset))
  }) + with(patient_data, {
    subset <- outlier_est[outlier]
    csubset <- outlier_est[outlier_det_2sigma]
    c(sum(!csubset), length(subset))
  })
  
  # estimation
  # fixed
  muhat_fixed <- muhat_fixed %>% arrange(time) %>% 
    left_join(unique(sim_data[sim_data$source == "Onsite", c("time", "y_true_fixed")]), by = "time")
  fwc_fixed <- muhat_fixed %>%
      summarise(
        est = tail(mu_fixed, 1),
        true = tail(y_true_fixed, 1),
        .groups = "drop"
      )
  # fixed + random
  tmp_data <- if (model_type == 4) oracle_data else sim_data
  insample_est <- c()
  insample_est["bias"] <- tmp_data %>% filter(source == "Onsite") %>% 
    select(c(patient, time, y_true_all, y_fitted_all)) %>% 
    mutate(bias = (y_fitted_all - y_true_all)) %>% pull(bias) %>% mean()
  
  insample_est["rbias"] <- tmp_data %>% filter(source == "Onsite") %>% 
    select(c(patient, time, y_true_all, y_fitted_all)) %>% 
    mutate(bias = (y_fitted_all - y_true_all)) %>% 
    group_by(patient) %>% 
    summarise(rbias = mean(bias)/mean(sim_data$y_true_all[sim_data$source == "Onsite"])) %>% 
    pull(rbias) %>% mean()
  
  insample_est["rmse"] <- tmp_data %>% filter(source == "Onsite") %>% 
    select(c(patient, time, y_true_all, y_fitted_all)) %>% 
    mutate(sq_error = (y_fitted_all - y_true_all)^2) %>% 
    group_by(time) %>% 
    summarise(rmse = sqrt(mean(sq_error))) %>% pull(rmse) %>% mean()
  
  insample_est["bias_fwc"] <- tmp_data %>% filter(source == "Onsite") %>% 
    group_by(patient) %>%
    arrange(patient, time) %>%
    summarise(
      true = tail(y_true_all, 1),
      fitted = tail(y_fitted_all, 1),
      .groups = "drop"
    ) %>% mutate(bias = fitted - true) %>% pull(bias) %>% mean()
  
  insample_est["rbias_fwc"] <- tmp_data %>% filter(source == "Onsite") %>% 
    group_by(patient) %>%
    arrange(patient, time) %>%
    summarise(
      true = tail(y_true_all, 1),
      fitted = tail(y_fitted_all, 1),
      .groups = "drop"
    ) %>% mutate(rbias = (fitted - true)/mean(true)) %>% pull(rbias) %>% mean()
  
  insample_est["rmse_fwc"] <- tmp_data %>% filter(source == "Onsite") %>% 
    group_by(patient) %>%
    arrange(patient, time) %>%
    summarise(
      true = tail(y_true_all, 1),
      fitted = tail(y_fitted_all, 1),
      .groups = "drop"
    ) %>% mutate(rmse = (fitted - true)^2) %>% pull(rmse) %>% mean() %>% sqrt()
  
  # prediction on a test set
  outsample_est <- c()
  outsample_est["bias"] <- test_data %>% 
    select(c(patient, time, y_true_all, y_pred)) %>% 
    mutate(bias = (y_pred - y_true_all)) %>% pull(bias) %>% mean()
  
  outsample_est["rbias"] <- test_data %>% 
    select(c(patient, time, y_true_all, y_pred)) %>% 
    mutate(bias = (y_pred - y_true_all)) %>% 
    group_by(patient) %>% 
    summarise(rbias = mean(bias)/mean(test_data$y_true_all)) %>% pull(rbias) %>% mean()
  
  outsample_est["rmse"] <- test_data %>% 
    select(c(patient, time, y_true_all, y_pred)) %>% 
    mutate(sq_error = (y_pred - y_true_all)^2) %>% 
    group_by(time) %>% 
    summarise(rmse = sqrt(mean(sq_error))) %>% pull(rmse) %>% mean()
  
  outsample_est["bias_fwc"] <- test_data %>% 
    group_by(patient) %>%
    arrange(patient, time) %>%
    summarise(
      true = tail(y_true_all, 1),
      fitted = tail(y_pred, 1),
      .groups = "drop"
    ) %>% mutate(bias = fitted - true) %>% pull(bias) %>% mean()
  
  outsample_est["rbias_fwc"] <- test_data %>% 
    group_by(patient) %>%
    arrange(patient, time) %>%
    summarise(
      true = tail(y_true_all, 1),
      fitted = tail(y_pred, 1),
      .groups = "drop"
    ) %>% mutate(rbias = (fitted - true)/mean(true)) %>% pull(rbias) %>% mean()
  
  outsample_est["rmse_fwc"] <- test_data %>% 
    group_by(patient) %>%
    arrange(patient, time) %>%
    summarise(
      true = tail(y_true_all, 1),
      fitted = tail(y_pred, 1),
      .groups = "drop"
    ) %>% mutate(rmse = (fitted - true)^2) %>% pull(rmse) %>% mean() %>% sqrt()

  ### save data
  save_var <- list(
    method      = method,
    convergence = convergence,
    data_fit    = res$data_fit,
    sim_data    = sim_data,
    test_data   = test_data,
    muhat_fixed = muhat_fixed,
    model_fit   = if (model_type %in% 1:5) res$model_fit else NULL,
    mu_fixed    = if (model_type == 6) res$mu_fixed else NULL,
    oracle_data = if (model_type == 4) oracle_data else NULL
  )
  saveRDS(save_var, file = paste0("./raw_results/raw_results_sc00_", rep, ".rds"))
  
  res_all <- tibble(
    method = method,
    replicate = rep,
    n_eff = convergence[1],
    Rhat = convergence[2], 
    converge_flag = convergence[3], 
    fpr = fpr,
    fnr = fnr,
    cfpr_3sigma = cfpr_3sigma, 
    cfpr_2sigma = cfpr_2sigma, 
    cfnr_3sigma = cfnr_3sigma,
    cfnr_2sigma = cfnr_2sigma,
    fpr_p_num = fpr_p[1], 
    fpr_p_denom = fpr_p[2],
    fnr_p_num = fnr_p[1], 
    fnr_p_denom = fnr_p[2],
    cfpr_p_num = cfpr_p[1], 
    cfpr_p_denom = cfpr_p[2], 
    error_p_num = error_p[1],
    error_p_denom = error_p[2],
    cerror_3sigma_p_num = cerror_3sigma_p[1],
    cerror_3sigma_p_denom = cerror_3sigma_p[2],
    cerror_2sigma_p_num = cerror_2sigma_p[1],
    cerror_2sigma_p_denom = cerror_2sigma_p[2],
    fwc_fixed_est = fwc_fixed$est,
    fwc_fixed_true = fwc_fixed$true,
    CI_length = CI_length,
    CI_cov = CI_cov,
    CI_length_fwc = CI_length_fwc,
    CI_cov_fwc = CI_cov_fwc,
    HPD_length = HPD_length,
    HPD_cov = HPD_cov,
    HPD_length_fwc = HPD_length_fwc,
    HPD_cov_fwc = HPD_cov_fwc,
    insp_bias = insample_est["bias"],
    insp_rbias = insample_est["rbias"],
    insp_rmse = insample_est["rmse"],
    insp_bias_fwc = insample_est["bias_fwc"],
    insp_rbias_fwc = insample_est["rbias_fwc"],
    insp_rmse_fwc = insample_est["rmse_fwc"],
    outsp_bias = outsample_est["bias"],
    outsp_rbias = outsample_est["rbias"],
    outsp_rmse = outsample_est["rmse"], 
    outsp_bias_fwc = outsample_est["bias_fwc"],
    outsp_rbias_fwc = outsample_est["rbias_fwc"],
    outsp_rmse_fwc = outsample_est["rmse_fwc"]
  )
  muhat_fixed_all <- bind_cols(
    tibble(method = method, replicate = rep),
    tibble(!!!setNames(muhat_fixed[["mu_fixed"]], paste0("time", muhat_fixed[["time"]])))
  )
  mutrue_fixed_all <- bind_cols(
    tibble(method = method, replicate = rep),
    tibble(!!!setNames(muhat_fixed[["y_true_fixed"]], paste0("time", muhat_fixed[["time"]])))
  )


  return(list(res_all = res_all, muhat_fixed_all = muhat_fixed_all, mutrue_fixed_all = mutrue_fixed_all))
}

### RUN code =======

output <- main_func(
  group = "trt",
  n_patients = all.config$n_patients[sc00],
  n_measure_o = 8,
  n_measure_r = 8,
  beta_delta = beta_delta_setting[all.config$beta_delta_idx[sc00], ],
  rho_source = all.config$rho_source[sc00],
  re_factor = all.config$re_factor[sc00],
  shift_rate = 0.3,
  remote_inflation = all.config$remote_inflation[sc00],
  outlier_rate1 = all.config$outlier_rate1[sc00],
  outlier_rate2 = all.config$outlier_rate2[sc00],
  outlier_type = all.config$outlier_type[sc00],
  uni_modal = all.config$uni_modal[sc00],
  high_freq = all.config$high_freq[sc00],
  DGM_type = all.config$DGM_type[sc00],
  model_type = all.config$model_type[sc00],
  gmm_ncls = all.config$gmm_ncls[sc00],
  gmm_type = all.config$gmm_type[sc00],
  error_type = all.config$error_type[sc00], 
  seed = seed00,
  rep = rep00
)

resultpath1 <- "./results/output_sc00.csv"
resultpath2 <- "./results/muhat_fixed_sc00.csv"
resultpath3 <- "./results/mutrue_fixed_sc00.csv"
### write file in a parallel way
lockfile <- "./results/lockfile.lock"
lock <- lock(lockfile, timeout = Inf)
if (!file.exists(resultpath1)) {
  write_excel_csv(output$res_all, file = resultpath1, append = T, col_names = TRUE)
} else {
  write_excel_csv(output$res_all, file = resultpath1, append = T, col_names = FALSE)
}

if (!file.exists(resultpath2)) {
  write_excel_csv(output$muhat_fixed_all, file = resultpath2, append = T, col_names = TRUE)
} else {
  write_excel_csv(output$muhat_fixed_all, file = resultpath2, append = T, col_names = FALSE)
}

if (!file.exists(resultpath3)) {
  write_excel_csv(output$mutrue_fixed_all, file = resultpath3, append = T, col_names = TRUE)
} else {
  write_excel_csv(output$mutrue_fixed_all, file = resultpath3, append = T, col_names = FALSE)
}
unlock(lock)
