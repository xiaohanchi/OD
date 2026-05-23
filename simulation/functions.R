
generate_data <- function(
    n_patients, n_measure = c("Onsite" = NA, "Remote" = NA),
    beta, beta_delta, sd_r, corr_r = 0.1, rho_source = 0.5, re_factor = 1, rho_ar1 = NULL, 
    sd_eps = NULL, shift_rate, remote_inflation, outlier_rate1, outlier_rate2,
    equal_remote_N = TRUE, outlier_type, uni_modal = TRUE, high_freq = FALSE,
    DGM_type, error_type, oracle = FALSE, seed = 123) {
  
  # sd_eps: sd for iid noise / AR error (u_{ij})
  # remote_inflation: inflation factor on sd for remote random effect (e.g., 0.2)
  # outlier_type: 1 for sporadic outliers; 2 for unimodal recurring outliers; 3 for mixed type
  # outlier_rate1: patient level (xx% of patients have outliers)
  # outlier_rate2: patient level (xx% of points for each patient are outliers)
  # uni_modal: used with outlier_type = 2 and 3
  # DGM_type: two types of data generating models
  # error_type: 1 for iid error; 2 for AR error
  
  high_freq_factor <- 8
  
  set.seed(seed)
  if(outlier_type == 3) {
    outlier_rate1 <- c(0.5, 0.2)
    outlier_rate2 <- c(0.2, 0.5)
  }
  
  remote_factor <- 2
  if (equal_remote_N) {
    remote_nmeasure_vec <- rep(n_measure["Remote"], n_patients) %>% unname()
  } else {
    remote_nmeasure_vec <- sapply(
      1:n_patients,
      function(rr) rpois(1, lambda = remote_factor * n_measure["Remote"])
    )
  }
  remote_nmeasure_vec <- pmax(remote_nmeasure_vec, 2) # at least two remote dta points
  
  if (outlier_type %in% c(1, 2)) {
    outlier_pts <- sample(n_patients, round(outlier_rate1 * n_patients), replace = FALSE)
    normal_pts <- setdiff(seq_len(n_patients), outlier_pts)
    if (outlier_type %in% c(2) & high_freq) {
      remote_nmeasure_vec[outlier_pts] <- sapply(
        1:length(outlier_pts),
        function(rr) rpois(1, lambda = high_freq_factor * n_measure["Remote"])
      )
      remote_nmeasure_vec[sample(normal_pts, round(length(normal_pts) * 0.5))] <- sapply(
        1:round(length(normal_pts) * 0.5),
        function(rr) rpois(1, lambda = high_freq_factor * n_measure["Remote"])
      )
    }
  } else {
    outlier_pts <- list(
      type1 = sample(n_patients, round(outlier_rate1[1] * n_patients), replace = FALSE),
      type2 = sample(n_patients, round(outlier_rate1[2] * n_patients), replace = FALSE)
      )
    normal_pts <- setdiff(seq_len(n_patients), unique(unlist(outlier_pts)))
    if (high_freq) {
      remote_nmeasure_vec[outlier_pts$type2] <- sapply(
        1:length(outlier_pts$type2),
        function(rr) rpois(1, lambda = high_freq_factor * n_measure["Remote"])
      )
      remote_nmeasure_vec[sample(normal_pts, round(length(normal_pts) * 0.5))] <- sapply(
        1:round(length(normal_pts) * 0.5),
        function(rr) rpois(1, lambda = high_freq_factor * n_measure["Remote"])
      )
    }
    outlier2_uni_pts <- sample(outlier_pts$type2, size = ceiling(0.5 * length(outlier_pts$type2)))  
  }

  
  # remote data: shift_rate*100% patients have shift in both random intercept and slope and (outlier_rate1, outlier_rate2)*100% outliers
  # time_final <- 10
  # time <- list(
  #   onsite = (0:(n_measure["Onsite"] - 1)) * ((time_final - 1) / (n_measure["Onsite"] - 1)),
  #   remote = lapply(
  #     remote_nmeasure_vec,
  #     function(nn) (0:(nn - 1)) * ((time_final - 1) / (nn - 1))
  #   )
  # )
  time <- list(
    onsite = 1:n_measure["Onsite"],
    remote = lapply(
      remote_nmeasure_vec,
      # function(nn) (0:(nn - 1)) * ((n_measure["Onsite"] - 1) / (nn - 1))
      function(nn) seq(1, n_measure["Onsite"], length.out = nn) + 0.05
    )
  )
  
  all_data <- rbind(
    data.frame(
      patient = rep(1:n_patients, each = n_measure["Onsite"]),
      time = rep(time$onsite, n_patients), source = "Onsite"
    ),
    data.frame(
      patient = rep(1:n_patients, times = remote_nmeasure_vec),
      time = unlist(time$remote), source = "Remote"
    )
  ) %>% mutate(
    r0 = NA, r1 = NA, y_true_fixed = NA, y_true_all = NA, y = NA, 
    eps_raw = NA, eps = NA, eps_flag_2sigma = FALSE, eps_flag_3sigma = FALSE, 
    outlier_undetect_2sigma = FALSE, outlier_undetect_3sigma = FALSE, 
    shift = FALSE, outlier1 = FALSE, outlier2 = FALSE
    )
  
  if (oracle) all_data <- all_data %>% arrange(patient, time)
  onsite_idx <- which(all_data$source == "Onsite")
  remote_idx <- which(all_data$source == "Remote")
  if (oracle) all_data$source <- "Onsite"

  # remote data + shift
  shift_pts <- sample(n_patients, round(shift_rate * n_patients), replace = FALSE)
  shift_idx <- which(all_data$patient %in% shift_pts & all_data$source == "Remote")
  if (!oracle) all_data$shift[shift_idx] <- TRUE
  
  # remote data + outliers
  if (!oracle) {
    if (outlier_type %in% c(1, 2)) {
      outlier_idx <- sapply(outlier_pts, function(id) {
        which(all_data$patient == id & all_data$source == "Remote") %>%
          sample(., round(outlier_rate2 * length(.)), replace = FALSE)
      }) %>% unlist()
      all_data[[paste0("outlier", outlier_type)]][outlier_idx] <- TRUE
    } else {
      for (ii in 1:2) {
        outlier_idx <- sapply(outlier_pts[[paste0("type", ii)]], function(id) {
          which(all_data$patient == id & all_data$source == "Remote") %>%
            sample(., round(outlier_rate2[ii] * length(.)), replace = FALSE)
        }) %>% unlist()
        all_data[[paste0("outlier", ii)]][outlier_idx] <- TRUE
      }
    }
  }
  
  
  
  Sigma_r <- matrix(
    c(sd_r[1]^2, (sd_r[1] * sd_r[2] * corr_r), (sd_r[1] * sd_r[2] * corr_r), sd_r[2]^2), 
    nrow = 2, byrow = TRUE
  )
  ns_basis <- ns(all_data$time, df = 3)
  
  set.seed(seed)
  for (ii in 1:n_patients) {
    r_onsite <- mvrnorm(n = 1, mu = rep(0, length(sd_r)), Sigma = Sigma_r)
    r_indep  <- mvrnorm(n = 1, mu = rep(0, length(sd_r)), Sigma = re_factor^2*Sigma_r)
    r_remote <- rho_source * r_onsite + sqrt(1 - rho_source^2) * r_indep
    # r_onsite <- mvrnorm(n = 1, mu = rep(0, length(sd_r)), Sigma = Sigma_r)
    # r_remote <- mvrnorm(n = 1, mu = rep(0, length(sd_r)), Sigma = Sigma_r)
    eps_raw <- if (error_type == 1) {
      list()
    } else {
      list(Onsite = rnorm(1, 0, sd_eps[1]), Remote = rnorm(1, 0, sd_eps[2]))
    }
    for (ss in c("Onsite", "Remote")) {
      if (oracle & ss ==  "Remote") next
      r_vals <- if (ss == "Onsite") r_onsite else r_remote
      
      tmp_idx <- which(all_data$patient == ii & all_data$source == ss)
      if (error_type == 1) {
        eps_raw[[ss]] <- rnorm(
          length(all_data$patient[tmp_idx]), 0, 
          sd_eps[match(ss, c("Onsite", "Remote"))]
          )
      } else if (error_type == 2) {
        noise_ar1 <- rnorm(
          length(all_data$patient[tmp_idx]), 0, 
          sd_eps[(1 + 1 * (ss == "Remote"))]
        )
        for(tt in 2:length(tmp_idx)) {
          if (ss == "Onsite") {
            eps_raw$Onsite[tt] <- rho_ar1[1] * eps_raw$Onsite[tt - 1] + noise_ar1[tt]
          } else {
            eps_raw$Remote[tt] <- rho_ar1[2] * 0.5 * (eps_raw$Remote[tt - 1] + eps_raw$Onsite[floor(all_data$time[tmp_idx][tt])]) + noise_ar1[tt]
          }
        }
      }
      eps <- eps_raw
      
      # random effects inflation for some patients
      if (all(all_data$shift[tmp_idx])) {
        r_vals <- mvrnorm(n = 1, mu = rep(0, length(sd_r)), Sigma = Sigma_r * (1 + remote_inflation)^2)
      }
      # + outlier for some points/patients
      # outlier_undetect <- rep(FALSE, length(tmp_idx))
      if (outlier_type == 1) {
        outlier1_flag <- all_data$outlier1[tmp_idx]
        if (sum(outlier1_flag) > 0) {
          eps[[ss]][outlier1_flag] <- rnorm(
            sum(outlier1_flag), 0, (10 * sd_eps[2])
          )
          # outlier_undetect[which(outlier1_flag & (abs(eps[[ss]]) < 3 * sd_eps[2]))] <- TRUE
        }
      } else if (outlier_type == 2) {
        outlier2_flag <- all_data$outlier2[tmp_idx]
        if (sum(outlier2_flag) > 0) {
          if (uni_modal) {
            outlier_mean <- rnorm(1, -60, 3)
            eps[[ss]][outlier2_flag] <- eps[[ss]][outlier2_flag] + outlier_mean
          } else {
            outlier_mean <- c(rnorm(1, -80, 3), rnorm(1, -40, 3))
            shuffled_idx <- sample(which(outlier2_flag))
            half <- length(shuffled_idx) %/% 2
            eps[[ss]][shuffled_idx[1:half]] <- eps[[ss]][shuffled_idx[1:half]] + outlier_mean[1]
            eps[[ss]][shuffled_idx[(half + 1):length(shuffled_idx)]] <- eps[[ss]][shuffled_idx[(half + 1):length(shuffled_idx)]] + outlier_mean[2]
          }
        }
        # outlier_undetect[which(outlier2_flag & (abs(eps[[ss]]) < 3 * sd_eps[2]))] <- TRUE
      } else if (outlier_type == 3) {
        outlier1_flag <- all_data$outlier1[tmp_idx]
        outlier2_flag <- all_data$outlier2[tmp_idx]
        if (sum(outlier1_flag) > 0) {
          eps[[ss]][(outlier1_flag & !outlier2_flag)] <- rnorm(
            sum(outlier1_flag & !outlier2_flag), 0, (10 * sd_eps[2])
          )
        }
        if (sum(outlier2_flag) > 0) {
          if (ii %in% outlier2_uni_pts) {
            outlier_mean <- rnorm(1, -60, 3)
            eps[[ss]][outlier2_flag] <- eps[[ss]][outlier2_flag] + outlier_mean 
          } else {
            outlier_mean <- c(rnorm(1, -80, 3), rnorm(1, -40, 3))
            shuffled_idx <- sample(which(outlier2_flag))
            half <- length(shuffled_idx) %/% 2
            eps[[ss]][shuffled_idx[1:half]] <- eps[[ss]][shuffled_idx[1:half]] + outlier_mean[1]
            eps[[ss]][shuffled_idx[(half + 1):length(shuffled_idx)]] <- eps[[ss]][shuffled_idx[(half + 1):length(shuffled_idx)]] + outlier_mean[2]
          }
        }
        # outlier_undetect[which((outlier1_flag | outlier2_flag) & (abs(eps[[ss]]) < 3 * sd_eps[2]))] <- TRUE
      }
      all_data$r0[tmp_idx] <- r_vals[1]
      all_data$r1[tmp_idx] <- r_vals[2]
      beta_val <- if (ss == "Onsite") beta else (beta + beta_delta)
      if (DGM_type == 1) {
        y_true_fixed <- if (length(tmp_idx) > 0) {
          sapply(all_data$time[tmp_idx], function(t) sum(beta_val * t^(0:(length(beta_val) - 1))))
        } else numeric(0)
      } else if (DGM_type == 2) {
        y_true_fixed <-  c(ns_basis[tmp_idx, ] %*% beta_val)
      }
      y_true_all <- if (length(tmp_idx) > 0) {
        y_true_fixed + sapply(all_data$time[tmp_idx], function(t) sum(r_vals * t^(0:(length(r_vals) - 1))))
      } else numeric(0)
      y_obs <- y_true_all + eps[[ss]]
      all_data$y_true_fixed[tmp_idx] <- y_true_fixed
      all_data$y_true_all[tmp_idx] <- y_true_all
      all_data$y[tmp_idx] <- y_obs
      all_data$eps_raw[tmp_idx] <- eps_raw[[ss]]
      all_data$eps[tmp_idx] <- eps[[ss]]
    }
  }
  
  # flag "outlier" normal points and "undetectable" outliers
  for (sigma_cutoff in c(3, 2)) {
    flag_col <- paste0("eps_flag_", sigma_cutoff, "sigma")
    undetect_col <- paste0("outlier_undetect_", sigma_cutoff, "sigma")
    for (ss in c("Onsite", "Remote")) {
      if (oracle & ss ==  "Remote") next
      
      if (error_type == 1) {
        sd_tmp <- sd_eps[match(ss, c("Onsite", "Remote"))]
      } else if (error_type == 2) {
        sd_tmp <- all_data %>%
          filter(source == ss, !outlier1, !outlier2) %>%
          summarise(sd_tmp = sd(eps_raw)) %>% pull()
      }
      
      idx_tmp1 <- with(
        all_data,
        which((abs(eps_raw) > sigma_cutoff * sd_tmp) & (source == ss))
      )
      all_data[[flag_col]][idx_tmp1] <- TRUE
      
      idx_tmp2 <- with(
        all_data,
        which((outlier1 | outlier2) & (abs(eps) < sigma_cutoff * sd_tmp))
      )
      all_data[[undetect_col]][idx_tmp2] <- TRUE
    }
  }
  
  if(oracle) all_data$source[remote_idx] <- "Remote"
  
  return(all_data)
}

twostage_lmm <- function(data_all, k = 5) {
  data_all$source <- factor(data_all$source)
  # fit onsite lmm
  fit.onsite <- lmer(y ~ ns(time, df = k) + (1 + ns(time, df = 3) | patient), 
                     data = (data_all %>% filter(source == "Onsite")))
  
  # fit remote lmm
  df.remote <- data_all %>% 
    filter(source == "Remote") %>%
    mutate(outlier_est = FALSE)
  fit.remote <- lmer(y ~ ns(time, df = k) + (1 + ns(time, df = 3) | patient), data = df.remote)
  df.remote$outlier_est[abs(residuals(fit.remote)) > 3 * sigma(fit.onsite)] <- TRUE
  
  data_fit <- data_all %>% left_join(df.remote[c("patient", "time", "y", "outlier_est")]) %>% 
    mutate(outlier_est = coalesce(outlier_est, FALSE))
  
  # final model fitting after excluding outliers
  fit.final <- lmer(y ~ ns(time, df = k) + (1 + ns(time, df = 3) | patient), 
                    data = filter(data_fit, !outlier_est))
  
  return(list(data_fit = data_fit, model_fit = fit.final))
}


run_lmm <- function(data_all, type, k = 5, error_type = 1) {
  # type 1: onsite data only; 2: pool onsite + remote data; 3: oracle model
  data_all$source <- factor(data_all$source)
  if (type == 1) {
    data_fit <- data_all %>% filter(source == "Onsite")
    # fit <- lmer(
    #   y ~ ns(time, df = k) + (1 + time | patient), data = data_fit
    # )
    fit <- lmer(
      y ~ ns(time, df = k) + (1 + ns(time, df = 3) | patient), data = data_fit
    )
  } else if (type == 2) {
    data_fit <- data_all
    # fit <- lmer(
    #   y ~ ns(time, df = k) + (1 + time | patient), data = data_fit
    # )
    fit <- lmer(
      y ~ ns(time, df = k) + (1 + ns(time, df = 3) | patient), data = data_fit
    )
  } else if (type == 3) {
    # oracle
    data_fit <- data_all
    if (error_type == 1) {
      fit <- lmer(
        y ~ poly(time, 2, raw = TRUE) + (1 + time || patient/source), data = data_fit
      )
    } else if (error_type == 2) {
      fit <- tryCatch({
        lme(
          y ~ poly(time, 2, raw = TRUE),
          random = list(patient = ~ pdDiag(~ 1 + time), source = ~ 1),
          correlation = corAR1(form = ~ time | patient/source),
          data = data_fit
        )
      }, error = function(e) {
        lme(
          y ~ poly(time, 2, raw = TRUE),
          random = list(patient = ~ 1, source = ~ 1),
          correlation = corAR1(form = ~ time | patient/source),
          data = data_fit
        )
      })
      message("Model used: ", deparse(fit$call$random))
    }
  }
  
  data_fit$outlier_est <- FALSE
  
  return(list(data_fit = data_fit, model_fit = fit))
}


run_robustlmm <- function(data_all) {
  data_all$source <- factor(data_all$source)
  fit <- rlmer(
    y ~ ns(time, df = 5) + (1 +  ns(time, df = 3) | patient),
    data = data_all # , weights = data_all$w
  )
  data_all$rweights <- getME(fit, "w_e")
  cutoff <- quantile(
    data_all$rweights[which(data_all$source == "Onsite")],
    0.05,
    na.rm = TRUE
  )
  
  data_all$outlier_est <- FALSE
  data_all$outlier_est[which(data_all$rweights < cutoff)] <- TRUE
  
  return(list(data_fit = data_all, model_fit = fit))
}

set_splines <- function(time, k) {
  knots <- quantile(unique(time), seq(0, 1, length = (k + 2))[-c(1, (k + 2))])
  Z_K <- (abs(outer(time, knots, "-")))^3
  OMEGA_all <- (abs(outer(knots, knots, "-")))^3
  svd.OMEGA_all <- svd(OMEGA_all)
  sqrt.OMEGA_all <- t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u) * sqrt(svd.OMEGA_all$d)))
  Z <- t(solve(sqrt.OMEGA_all, t(Z_K)))
  
  return(list(
    Z = Z,
    knots = knots,
    sqrt.OMEGA_all = sqrt.OMEGA_all
  ))
}

predict_splines <- function(new_time, spline_info) {
  knots <- spline_info$knots
  sqrt.OMEGA_all <- spline_info$sqrt.OMEGA_all
  
  Z_K_new <- (abs(outer(new_time, knots, "-")))^3
  Z_new <- t(solve(sqrt.OMEGA_all, t(Z_K_new)))
  
  return(Z_new)
}

run_gmm <- function(data_all, test_data, n_cls = 5, var_inflate = 5,
                    a0 = 10, a_slab, b_slab, type, k1 = 5, k2 = 3) {
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  # var_inflate: variance inflation factor for the outliers
  df.fit <- data_all %>% rename(id = patient, y_obs = y, group = source)
  
  df.fit$newid <- as.integer(factor(df.fit$id))
  
  B_fit1 <- set_splines(time = df.fit$time, k = k1)
  B_est1 <- predict_splines(new_time = test_data$time, B_fit1)
  B_fit2 <- set_splines(time = df.fit$time, k = k2)
  B_est2 <- predict_splines(new_time = test_data$time, B_fit2)
  
  if (type %in% c(1:2)) {
    # type = 1,2: discarded
  } else if (type %in% c(3:5)) {
    mcmc_data <- list(
      n_cls = n_cls, var_inflate = var_inflate, a0 = a0,
      y = df.fit$y_obs, time = df.fit$time,
      group = ifelse(df.fit$group == "Onsite", 1, 0),
      N = nrow(df.fit), n_patients = length(unique(df.fit$newid)),
      patient_id = df.fit$newid,
      n_ref = sum(df.fit$group == "Onsite"),
      n_rem = sum(df.fit$group != "Onsite"),
      idx_ref = which(df.fit$group == "Onsite"),
      idx_rem = which(df.fit$group != "Onsite"),
      N_test = nrow(test_data), 
      n_patients_test = length(unique(test_data$patient)),
      time_test = test_data$time, 
      patient_id_test = test_data$patient,
      B1 = B_fit1$Z, B_test1 = B_est1, K1 = k1,
      B2 = B_fit2$Z, B_test2 = B_est2, K2 = k2
    )
  }
  
  paras <- c("beta0", "beta1", "b_spline", "sigma_ref", "sigma", "z", "mu_fixed", "mu_y", "y_pred")
  
  # initial values for "parameters" block
  # init_val <- list(sigma_ref = 1.5, sigma_rem_d = 1.5)
  init_fun <- function() {
    list(sigma_ref = 4, sigma_rem_d = 2)
  }
  
  stan_fit <- stan(
    model_code = case_when(
      type == 3 ~ sGMM.pool,
      type == 4 ~ sGMM.pool.heteroRE,
      type == 5 ~ (str_replace_all(sGMM.heteroRE, "a_slab", a_slab) %>% 
        str_replace_all("b_slab", b_slab))
    ),
    data = mcmc_data,
    chains = 3,
    iter = 3500,
    warmup = 1000,
    thin = 2,
    seed = 123,
    init = init_fun,
    pars = paras,
    include = TRUE,
    control = list(adapt_delta = 0.85, max_treedepth = 10),
    verbose = FALSE, refresh = 200
  )
  
  time_onsite <- unique(df.fit$time[df.fit$group == "Onsite"])
  B_onsite <- predict_splines(new_time = time_onsite, B_fit1)
  if (type %in% c(1:2)) {

  } else if (type %in% c(3:5)) {
    z.mean <- apply(as.matrix(As.mcmc.list(stan_fit, "z")), 2, function(x) {
      tab <- table(factor(x, levels = 1:n_cls))
      tab / length(x)
    })
    mu.all <- As.mcmc.list(stan_fit, "mu_y") %>% as.matrix() %>% colMeans()
    y.pred <- As.mcmc.list(stan_fit, "y_pred") %>% as.matrix() %>% colMeans()
    df.fit$rweights <- t(z.mean)
    
    beta <- As.mcmc.list(stan_fit, c("beta0", "beta1")) %>% as.matrix()
    b_spline <- As.mcmc.list(stan_fit, c("b_spline")) %>% as.matrix()
    mu.fixed <- matrix(NA, nrow = nrow(beta), ncol = length(time_onsite))
    for (i in 1:nrow(beta)) {
      mu.fixed[i, ] <- sapply(time_onsite, function(t) sum(beta[i, ] * t^(0:(length(beta[i, ]) - 1)))) + B_onsite %*% b_spline[i, ]
    }
    # mu.fixed <- colMeans(mu.fixed)
  }
  df.fit$outlier_est <- FALSE
  if (is.null(dim(df.fit$rweights))) {
    outlier <- which(df.fit$rweights > 0.5)
  } else {
    outlier <- which(apply(df.fit$rweights, 1, which.max) != 1)
  }
  df.fit$outlier_est[outlier] <- TRUE
  
  return(list(data_fit = df.fit, stan_fit = stan_fit, 
              mu_fixed = mu.fixed, mu_all = mu.all, y_pred = y.pred))
}
