rm(list = ls())

library(pacman)
p_load(
  brms, coda, dplyr, readr, tidyr, tibble, lme4, mfp, nlme, robustlmm, 
  rjags, runjags, rstan, splines, stringr, mmrm, merTools, emmeans, filelock
)

load("data_perturbed2.RData")
source("stan_more.R")
### LMM ===========================
run_lmm <- function(data_all, trt.group, type, spline.type = 1) {
  df <- data.frame(
    id = data_all$USUBJID, time = data_all$time,
    y_obs = data_all$diffw, baseline = data_all$WEIGHTBL, 
    treatment = data_all$TRT01P, group = data_all$source
  )
  # type 1: onsite data only; 2: pool onsite + remote data; 3: remote data only
  df$group <- factor(df$group)
  
  if (trt.group == 1) {
    df.fit <- df[which(df$treatment == "TZP MTD"), ]
  } else if (trt.group == 0) {
    df.fit <- df[which(df$treatment != "TZP MTD"), ]
  }
  
  if (type == 1) {
    df.fit <- df.fit %>% filter(group == "Onsite")
  } else if (type == 2) {
    df.fit <- df.fit
  } else if (type == 3) {
    df.fit <- df.fit %>% filter(group == "Remote")
  }
  
  if (spline.type == 1) {
    fit <- lmer(y_obs ~ baseline + ns(time, df = 4) + (1 + time | id), data = df.fit)
  } else {
    fit <- lmer(y_obs ~ baseline + ns(time, df = 4) + (1 | id) + (0 + ns(time, df = 2) | id), data = df.fit)
  }
  df.fit <- df.fit %>% mutate(y_fitted = NA_real_)
  df.fit$y_fitted[!is.na(df.fit$y_obs)] <- fitted(fit)

  return(list(df.fit = df.fit, model_fit = fit))
}

### Bayesian LMM ===========================
run_blmm <- function(data_all, trt.group, type, spline.type = 1, seed = 123) {
  
  df <- data.frame(
    id = data_all$USUBJID, time = data_all$time,
    y_obs = data_all$diffw, baseline = data_all$WEIGHTBL,
    treatment = data_all$TRT01P, group = data_all$source
  )
  df$group <- factor(df$group)
  
  if (trt.group == 1) {
    df.fit <- df[which(df$treatment == "TZP MTD"), ]
  } else if (trt.group == 0) {
    df.fit <- df[which(df$treatment != "TZP MTD"), ]
  }
  
  if (type == 1) {
    df.fit <- df.fit %>% filter(group == "Onsite")
  } else if (type == 2) {
    df.fit <- df.fit
  } else if (type == 3) {
    df.fit <- df.fit %>% filter(group == "Remote")
  }
  
  model_prior <- c(
    prior(normal(0, 50), class = Intercept),
    prior(normal(0, 50), class = b),        # fixed effects (baseline + spline)
    prior(student_t(3, 0, 5), class = sd),      # random effect SDs
    prior(student_t(3, 0, 5), class = sigma),   # residual SD
    prior(lkj(1), class = cor)              # correlation of random effects: flat prior
  )
  
  formula <- if (spline.type == 1) {
    bf(y_obs ~ baseline + ns(time, df = 4) + (1 + time | id))
  } else {
    bf(y_obs ~ baseline + ns(time, df = 4) + (1  | id) + (0 + ns(time, df = 2) | id))
  }
  
  fit <- brm(
    formula,
    data = df.fit,
    family = gaussian(),
    prior = model_prior,
    iter = 5000,
    warmup = 2000,
    chains = 3,
    cores = 3,
    thin = 2,
    seed = seed,
    control = list(adapt_delta = 0.85)
  )
  
  # fitted() "Estimate": posterior mean
  df.fit <- df.fit %>% mutate(y_fitted = NA_real_)
  df.fit$y_fitted[!is.na(df.fit$y_obs)] <- fitted(fit)[, "Estimate"]
  
  return(list(df.fit = df.fit, model_fit = fit))
}



### Mixture Model ============================

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

run.gmm <- function(data_all, trt.group, n_cls = 5, var_inflate = 10, 
                    a0 = 10, type, k1 = 5, k2 = 5, 
                    sigma_prior = "student_t(1, 0, 5)", 
                    lambda_prior = "student_t(1, 0, 1)",
                    tau_prior = "student_t(1, 0, 1)", 
                    a_slab = "0.1", b_slab = "0.1", seed = 123) {
  # type = 1: no onsite data (for outlier detection); 
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  # var_inflate: variance inflation factor for the outliers
  df <- data.frame(
    id = data_all$USUBJID, time = data_all$time,
    y_obs = data_all$diffw, baseline = data_all$WEIGHTBL, 
    treatment = data_all$TRT01P, group = data_all$source
  )

  if (trt.group == 1) {
    xbar_bl <- mean(data_baseline$WEIGHTBL[data_baseline$TRT01P == "TZP MTD"])
    df.fit <- df[which(df$treatment == "TZP MTD"), ]
    data_validation <- data_baseline_onsite %>% filter(TRT01P == "TZP MTD")
  } else if (trt.group == 0) {
    xbar_bl <- mean(data_baseline$WEIGHTBL[data_baseline$TRT01P != "TZP MTD"])
    df.fit <- df[which(df$treatment != "TZP MTD"), ]
    data_validation <- data_baseline_onsite %>% filter(TRT01P != "TZP MTD")
  }

  if (type %in% c(1:2)) {
    df.fit <- df.fit %>% filter(group == "Onsite")
  }
  df.fit$newid <- as.integer(factor(df.fit$id))
  
  df.fit <- df.fit %>% arrange(newid, group, time) %>% 
    group_by(newid, group) %>%
    mutate(ar_idx = row_number())
  
  time_est <- c(0:12)
  B_fit1 <- set_splines(time = df.fit$time, k = k1)
  B_est1 <- predict_splines(new_time = time_est, B_fit1)
  B_fit2 <- set_splines(time = df.fit$time, k = k2)
  B_est2 <- predict_splines(new_time = time_est, B_fit2)
  if (type %in% c(1)) {
    # Onsite GMM
    mcmc_data <- list(
      y = df.fit$y_obs, x_bl = df.fit$baseline, time = df.fit$time,
      N = nrow(df.fit), N_test = nrow(data_validation), 
      n_patients = length(unique(df.fit$newid)),
      patient_id = df.fit$newid, x_bl_test = data_validation$WEIGHTBL, 
      B = B_fit1$Z, B_test = B_est1[c(1, 13), ], K = k1
    )
    paras <- c("beta0", "beta1", "beta_bl", "b_spline", "r0", "r1","sigma_ref", "sigma_r0", "sigma_r1", "mu_fixed_test", "y_pred_t12")
    init_fun <- function() {
      list(sigma_ref = 2)
    }
  } else if (type %in% c(2)) {
    # Onsite GMM v2
    mcmc_data <- list(
      N_obs = sum(!is.na(df.fit$y_obs)), N_mis = sum(is.na(df.fit$y_obs)),
      idx_obs = which(!is.na(df.fit$y_obs)), idx_mis = which(is.na(df.fit$y_obs)),
      y_obs = df.fit$y_obs[!is.na(df.fit$y_obs)],
      x_bl = df.fit$baseline, time = df.fit$time,
      N = nrow(df.fit), N_test = nrow(data_validation),
      n_patients = length(unique(df.fit$newid)),
      patient_id = df.fit$newid, x_bl_test = data_validation$WEIGHTBL,
      B1 = B_fit1$Z, B_test1 = B_est1[c(1, 13), ], K1 = k1,
      B2 = B_fit2$Z, B_test2 = B_est2[c(1, 13), ], K2 = k2
    )
    paras <- c("beta0", "beta1", "beta_bl", "b_spline", "r0", "r1","sigma_ref", "sigma_r0", "sigma_r1", "mu_fixed_test", "y_pred_t12")
    init_fun <- function() {
      list(sigma_ref = 2)
    }
  } else if (type %in% c(3)) {
    # Pooled sGMM
    mcmc_data <- list(
      n_cls = n_cls, var_inflate = var_inflate, a0 = a0, 
      y = df.fit$y_obs, time = df.fit$time,
      x_bl = df.fit$baseline, time = df.fit$time,
      N = nrow(df.fit), N_test = nrow(data_validation), 
      n_patients = length(unique(df.fit$newid)),
      patient_id = df.fit$newid, 
      x_bl_test = data_validation$WEIGHTBL, 
      n_ref = sum(df.fit$group == "Onsite"),
      n_rem = sum(df.fit$group != "Onsite"),
      idx_ref = which(df.fit$group == "Onsite"), 
      idx_rem = which(df.fit$group != "Onsite"), 
      B = B_fit1$Z, B_test = B_est1[c(1, 13), ], K = k1
    )
    paras <- c("beta0", "beta1",  "beta_bl", "b_spline", "r0", "r1","sigma_ref", "sigma", "sigma_r0", "sigma_r1", "z", "mu_fixed_test", "y_pred_t12")
    init_fun <- function() {
      list(sigma_ref = 2, sigma_rem = 2)
    }
    
  } else if (type %in% c(4:5)) {
    # Pooled sGMM v2
    mcmc_data <- list(
      n_cls = n_cls, var_inflate = var_inflate, a0 = a0,
      N_obs = sum(!is.na(df.fit$y_obs)), N_mis = sum(is.na(df.fit$y_obs)),
      idx_obs = which(!is.na(df.fit$y_obs)), idx_mis = which(is.na(df.fit$y_obs)),
      y_obs = df.fit$y_obs[!is.na(df.fit$y_obs)], time = df.fit$time,
      x_bl = df.fit$baseline, time = df.fit$time,
      N = nrow(df.fit), N_test = nrow(data_validation),
      n_patients = length(unique(df.fit$newid)),
      patient_id = df.fit$newid,
      x_bl_test = data_validation$WEIGHTBL,
      n_ref = sum(df.fit$group == "Onsite"),
      n_rem = sum(df.fit$group != "Onsite"),
      idx_ref = which(df.fit$group == "Onsite"),
      idx_rem = which(df.fit$group != "Onsite"),
      B1 = B_fit1$Z, B_test1 = B_est1[c(1, 13), ], K1 = k1,
      B2 = B_fit2$Z, B_test2 = B_est2[c(1, 13), ], K2 = k2
    )
    paras <- c("beta0", "beta1",  "beta_bl", "b_spline", "r0", "r1","sigma_ref", "sigma", "sigma_r0", "sigma_r1", "z", "mu_fixed_test", "y_pred_t12")
    init_fun <- function() {
      list(sigma_ref = 2, sigma_rem = 2)
    }
    
  } else if (type %in% c(6:12)) {
    mcmc_data <- list(
      n_cls = n_cls, var_inflate = var_inflate, a0 = a0,
      N_obs = sum(!is.na(df.fit$y_obs)), N_mis = sum(is.na(df.fit$y_obs)),
      idx_obs = which(!is.na(df.fit$y_obs)), idx_mis = which(is.na(df.fit$y_obs)),
      y_obs = df.fit$y_obs[!is.na(df.fit$y_obs)], time = df.fit$time,
      x_bl = df.fit$baseline, time = df.fit$time,
      N = nrow(df.fit), N_test = nrow(data_validation),
      n_patients = length(unique(df.fit$newid)),
      patient_id = df.fit$newid,
      x_bl_test = data_validation$WEIGHTBL,
      n_ref = sum(df.fit$group == "Onsite"),
      n_rem = sum(df.fit$group != "Onsite"),
      idx_ref = which(df.fit$group == "Onsite"),
      idx_rem = which(df.fit$group != "Onsite"),
      B1 = B_fit1$Z, B_test1 = B_est1[c(1, 13), ], K1 = k1,
      B2 = B_fit2$Z, B_test2 = B_est2[c(1, 13), ], K2 = k2
    )
    paras <- c("beta0", "beta1",  "beta_bl", "b_spline", "r0", "r1","sigma_ref", "sigma", "sigma_r0", "sigma_r1", "z", "mu_fixed_test", "y_pred_t12")
    init_fun <- function() {
      list(sigma_ref = 2, sigma_rem = 2)
    }
    
  }
  
  
  stan_fit <- stan(
    model_code = case_when(
      type == 1 ~ sGMM.onsite,
      type == 2 ~ sGMM.onsite.v2,
      type == 3 ~ sGMM.pool,
      type == 4 ~ sGMM.pool.v2,
      type == 5 ~ sGMM.pool.v2.hetero,
      type == 6 ~ (str_replace_all(sGMM.transfer6, "sigma_prior", sigma_prior) %>%
                     str_replace_all("lambda_prior", lambda_prior) %>% 
                     str_replace_all("tau_prior", tau_prior)),
      type == 7 ~ (str_replace_all(sGMM.transfer6.semisep, "lambda_prior", lambda_prior) %>% 
                      str_replace_all("tau_prior", tau_prior)),
      type == 8 ~ (str_replace_all(sGMM.transfer6.semisep2, "lambda_prior", lambda_prior) %>% 
                     str_replace_all("tau_prior", tau_prior)),
      type == 9 ~ (str_replace_all(sGMM.transfer6.semisep3, "lambda_prior", lambda_prior) %>% 
                     str_replace_all("tau_prior", tau_prior)),
      type == 10 ~ (str_replace_all(sGMM.transfer6.semisep4, "lambda_prior", lambda_prior) %>% 
                     str_replace_all("tau_prior", tau_prior)),
      type == 11 ~ (str_replace_all(sGMM.transfer6.semisep5.nocor, "lambda_prior", lambda_prior) %>% 
                      str_replace_all("tau_prior", tau_prior) %>% 
                      str_replace_all("a_slab", a_slab) %>% str_replace_all("b_slab", b_slab)),
      type == 12 ~ (str_replace_all(sGMM.transfer6.semisep5, "lambda_prior", lambda_prior) %>% 
                      str_replace_all("tau_prior", tau_prior) %>% 
                      str_replace_all("a_slab", a_slab) %>% str_replace_all("b_slab", b_slab)) 
      # type == 8 ~ (str_replace_all(sGMM.transfer6.heteroRE, "sigma_prior", sigma_prior) %>% 
      #                str_replace_all("lambda_prior", lambda_prior) %>% 
      #                 str_replace_all("tau_prior2", tau_prior) %>% str_replace_all("tau_prior", tau_prior)),
      # type == 9 ~ (str_replace_all(sGMM.transfer6.heteroRE2, "sigma_prior", sigma_prior) %>% 
      #                 str_replace_all("lambda_prior", lambda_prior) %>% 
      #                 str_replace_all("tau_prior2", tau_prior) %>% 
      #                 str_replace_all("tau_prior", tau_prior)),
      # type == 10 ~ (str_replace_all(sGMM.transfer6.heteroRE2, "sigma_prior", sigma_prior) %>% 
      #                 str_replace_all("lambda_prior", lambda_prior) %>% 
      #                 str_replace_all("tau_prior2", "student_t(1, 0, 0.1)") %>% 
      #                 str_replace_all("tau_prior", tau_prior))
    ),
    data = mcmc_data,
    chains = 3,
    iter = 6500, # 4000,
    warmup = 1500, # 1000,
    thin = 4,
    seed = seed,
    init = init_fun,
    pars = paras, 
    include = TRUE,
    control = list(adapt_delta = 0.85, max_treedepth = 10),
    verbose = FALSE , refresh = 200
  )

  if (!(type %in% c(1:2))) {
    z.mean <- apply(as.matrix(As.mcmc.list(stan_fit, "z")), 2, function(x) {
      tab <- table(factor(x, levels = 1:n_cls))
      tab / length(x)
    })
    df.fit$rweights <- t(z.mean)
  }
  
  beta <- As.mcmc.list(stan_fit, c("beta0", "beta1", "beta_bl")) %>% as.matrix()
  b_spline <- As.mcmc.list(stan_fit, c("b_spline")) %>% as.matrix()
  mu.fixed <- matrix(NA, nrow = nrow(beta), ncol = length(time_est))
  for (i in 1:nrow(beta)) {
    mu.fixed[i, ] <- sapply(time_est, function(t) sum(beta[i, 1:2] * t^(0:(length(beta[i, 1:2]) - 1)))) + B_est1 %*% b_spline[i, ] + beta[i, 3] * xbar_bl
  }

  return(list(df.fit = df.fit, stan_fit = stan_fit, mu.fixed = mu.fixed))
}

main.func <- function(data_all, trt.group, n_cls, var_inflate, 
                      a0, model.type, lmm.type, gmm.type, k1, k2, data.type, 
                      sigma_prior, lambda_prior, tau_prior, a_slab = 0.1, b_slab = 0.1, 
                      mcmc_seed = 123) {
  # data.type: 1 = full data; 2 = reduced data
  set.seed(233)
  
  if(data.type == 1) {
    # keep all visits; add NA's for any missing visit
    nom_time <- data_all %>%
      filter(source == "Onsite") %>%
      group_by(VISIT) %>%
      summarise(nom_time = round(median(time), 3), .groups = "drop") %>%
      arrange(nom_time) %>%
      deframe()
    data_run <- data_all %>%
      filter(source == "Onsite") %>%
      complete(USUBJID, VISIT = names(nom_time)) %>%
      group_by(USUBJID) %>%
      fill(SSTUDYID, TRT01P, WEIGHTBL, source, .direction = "downup") %>%
      ungroup() %>%
      mutate(time = coalesce(time, nom_time[VISIT])) %>%
      bind_rows(data_all %>% filter(source == "Remote"))
  } else if (data.type == 2) {
    # week 0, 12, 24, 36, 52, add NA's for missing target visits
    target_visits <- c("VISIT2", "VISIT5", "VISIT8", "VISIT9", "VISIT11")
    visit_times <- data_all %>% 
      filter(VISIT %in% c("VISIT2","VISIT5","VISIT8","VISIT9","VISIT11")) %>%
      group_by(VISIT, VISITNUM) %>%
      summarise(median_time = median(time), .groups = "drop") %>%
      arrange(VISITNUM) %>% pull(median_time)
    nom_time     <- setNames(round(visit_times, 3), target_visits)
    data_remote <- data_all %>% filter(source == "Remote")
    data_run <- data_all %>%
      filter(source == "Onsite") %>%
      filter(VISIT %in% target_visits) %>%
      complete(USUBJID, VISIT = target_visits) %>%
      group_by(USUBJID) %>%
      fill(SSTUDYID, TRT01P, WEIGHTBL, source, .direction = "downup") %>%
      ungroup() %>%
      mutate(
        time     = coalesce(time,     nom_time[VISIT]),
      ) %>%
      bind_rows(data_remote)
    
    # data_counts <- data_run %>% filter(source == "Onsite" & TRT01P != "Placebo") %>% 
    #   group_by(USUBJID) %>%
    #   summarise(
    #     n_visits = sum(!is.na(diffw)),
    #     .groups = "drop"
    #   )
    # table(data_counts$n_visits)
    
  } else if (data.type == 3) {
    # first, last, and week 24 (visit8) measurements 
    data_remote <- data_all %>% filter(source == "Remote") 
    data_run <- data_all %>%
      filter(source == "Onsite") %>%
      group_by(USUBJID) %>%
      arrange(USUBJID, ady, time) %>%
      mutate(
        time_max = max(ady),
        time_mid = time_max / 2,
        is_mid = row_number() == which.min(abs(ady - time_mid))
      ) %>%
      filter(row_number() == 1 | row_number() == n() | is_mid) %>%
      select(-time_max, -time_mid, -is_mid) %>%
      ungroup() %>%
      bind_rows(data_remote)
  }
  
  if (trt.group == 1) {
    xbar_bl <- mean(data_baseline$WEIGHTBL[data_baseline$TRT01P == "TZP MTD"])
    data_validation <- data_baseline_onsite %>% filter(TRT01P == "TZP MTD")
    y_obs <- data_onsite_only %>% filter(VISITNUM == 11 & TRT01P == "TZP MTD") %>% pull(diffw)
  } else if (trt.group == 0) {
    xbar_bl <- mean(data_baseline$WEIGHTBL[data_baseline$TRT01P != "TZP MTD"])
    data_validation <- data_baseline_onsite %>% filter(TRT01P != "TZP MTD")
    y_obs <- data_onsite_only %>% filter(VISITNUM == 11 & TRT01P != "TZP MTD") %>% pull(diffw)
  }
  
  if (model.type == 1) {
    # onsite lmm
    output <- run_lmm(data_all = data_run, trt.group = trt.group, type = 1, spline.type = lmm.type)
  } else if (model.type == 2) {
    # pool
    output <- run_lmm(data_all = data_run, trt.group = trt.group, type = 2, spline.type = lmm.type)
  } else if (model.type == 3) {
    # remote
    output <- run_lmm(data_all = data_run, trt.group = trt.group, type = 3, spline.type = lmm.type)
  } else if (model.type == 4) {
    # onsite Bayesian lmm
    output <- run_blmm(
      data_all = data_run, trt.group = trt.group, type = 1, spline.type = lmm.type, seed = mcmc_seed
      )
  } else if (model.type == 5) {
    # pool Bayesian
    output <- run_blmm(
      data_all = data_run, trt.group = trt.group, type = 2, spline.type = lmm.type, seed = mcmc_seed
    )
  } else if (model.type == 6) {
    # remote Bayesian
    output <- run_blmm(
      data_all = data_run, trt.group = trt.group, type = 3, spline.type = lmm.type, seed = mcmc_seed
    )
  } else if (model.type == 7) {
    output <- run.gmm(
      data_all = data_run, trt.group = trt.group, n_cls = n_cls, 
      var_inflate = var_inflate, a0 = a0, type = gmm.type, k1 = k1, k2 = k2, 
      sigma_prior = sigma_prior, lambda_prior = lambda_prior, tau_prior = tau_prior,
      a_slab = a_slab, b_slab = b_slab, seed = mcmc_seed
    )
  }
  
  if(model.type %in% c(1:3)) {
    # fixed effects
    CI_res_fwc <- emmeans(output$model_fit, ~ time, 
                          at = list(baseline = xbar_bl, time = c(0:12)),
                          params = "k", data = output$df.fit, level = 0.95, df = Inf) %>% 
      confint()
    est_fwc <- CI_res_fwc$emmean
    se_fwc <- CI_res_fwc$SE
    CI_l <- CI_res_fwc$asymp.LCL
    CI_u <- CI_res_fwc$asymp.UCL
    eps <- sigma(output$model_fit)
    convergence <- ""
    
    emm_res <- emmeans(output$model_fit, ~ time, 
                       at = list(baseline = xbar_bl, time = seq(0, 12, 1)),
                       params = "k", 
                       data = output$df.fit, 
                       level = 0.95, 
                       df = Inf)
    weights <- c(0.5, rep(1, length(seq(0, 12, 1)) - 2), 0.5) 
    auc_result <- contrast(emm_res, 
                           method = list("AUC" = weights),
                           adjust = "none") %>% confint()
    est_auc <- auc_result$estimate
    se_auc <- auc_result$SE
    ci_auc <- c(auc_result$asymp.LCL, auc_result$asymp.UCL)
    
    r0_sd <- attr(summary(output$model_fit)$varcor$id, "stddev")[1] %>% unname()
    r1_sd <- NA
    
  } else if (model.type %in% c(4:6)) {
    # fixed effects
    CI_res_fwc <- emmeans(output$model_fit, ~ time, 
                          at = list(baseline = xbar_bl, time = c(0:12)),
                          data = output$df.fit, level = 0.95, df = Inf) %>% 
      confint()
    est_fwc <- CI_res_fwc$emmean
    se_fwc <- (CI_res_fwc$upper.HPD - CI_res_fwc$lower.HPD) / (2 * qnorm(0.975))
    CI_l <- CI_res_fwc$lower.HPD
    CI_u <- CI_res_fwc$upper.HPD
    eps <- posterior_summary(output$model_fit, variable = "sigma")[, "Estimate"]
    rhat_vals  <- rhat(output$model_fit)
    convergence <- ifelse(all(rhat_vals < 1.02, na.rm = TRUE), "Yes", "No")
    
    emm_res <- emmeans(output$model_fit, ~ time, 
                       at = list(baseline = xbar_bl, time = seq(0, 12, 1)),
                       data = output$df.fit, level = 0.95, df = Inf)
    weights <- c(0.5, rep(1, length(seq(0, 12, 1)) - 2), 0.5) 
    auc_result <- contrast(emm_res, 
                           method = list("AUC" = weights),
                           adjust = "none") %>% confint()
    est_auc <- auc_result$estimate
    se_auc <- (auc_result$upper.HPD - auc_result$lower.HPD) / (2 * qnorm(0.975))
    ci_auc <- c(auc_result$lower.HPD, auc_result$upper.HPD)
    
    r0_sd <- summary(output$model_fit)$random$id[1, "Estimate"]
    r1_sd <- NA
    
  } else if (model.type == 7) {
    # fixed
    est_fwc <- colMeans(output$mu.fixed)
    se_fwc <- apply(output$mu.fixed, 2, sd)
    quantiles <- apply(output$mu.fixed, 2, function(x) {
      quantile(x, probs = c(0.025, 0.975))
    })
    CI_l <- quantiles[1, ]
    CI_u <- quantiles[2, ]
    eps <- NA
    rhat <- tryCatch(
      summary(output$stan_fit)$summary[c("beta0", "beta1", "sigma_ref", "sigma[1]"), "Rhat"],
      error = function(e) summary(output$stan_fit)$summary[c("beta0", "beta1", "sigma_ref"), "Rhat"]
      ) %>% round(2) %>% max()
    convergence <- ifelse(rhat > 1.03, "No", "Yes")
    
    weights <- c(0.5, rep(1, length(seq(0, 12, 1)) - 2), 0.5) 
    auc_samples <- apply(output$mu.fixed, 1, function(mu_sample) {
      sum(weights * mu_sample)
    })
    est_auc <- mean(auc_samples)
    se_auc <- sd(auc_samples)
    ci_auc <- quantile(auc_samples, probs = c(0.025, 0.975))
    
    r0_sd <- summary(output$stan_fit)$summary[c("sigma_r0"), "mean"]
    r1_sd <- summary(output$stan_fit)$summary[c("sigma_r1"), "mean"]
    
    # #random
    # samples_mat <- As.mcmc.list(output$stan_fit, c("y_pred_t12")) %>% as.matrix()
    # get_re_summary <- function(mat, par_name, y_obs) {
    #   cols <- mat[, grep(paste0("^", par_name, "\\["), colnames(mat))]
    #   id_means <- colMeans(cols)
    #   id_sds <- apply(cols, 2, sd)
    #   pred_CI <- apply(cols, 2, quantile, probs = c(0.025, 0.975))
    #   tibble(
    #     mean_abs = mean(id_means), 
    #     median_sd = median(id_sds), 
    #     CI_len = mean(pred_CI[2, ] - pred_CI[1, ]), 
    #     coverage = mean(y_obs <= pred_CI[2, ] & y_obs >= pred_CI[1, ]) * 100
    #   )
    # }
    # pred_full <- get_re_summary(samples_mat, "y_pred_t12", y_obs = y_obs)
    # 
    # #pred
    # samples_mat <- As.mcmc.list(output$stan_fit, c("mu_fixed_test")) %>% as.matrix()
    # get_re_summary <- function(mat, par_name) {
    #   cols <- mat[, grep(paste0("^", par_name, "\\["), colnames(mat))]
    #   id_means <- colMeans(cols)
    #   id_sds <- apply(cols, 2, sd)
    #   tibble(
    #     pred_bias = mean(id_means - y_obs), 
    #     pred_rmse = sqrt(mean((id_means - y_obs)^2))
    #   )
    # }
    # pred_res <- get_re_summary(samples_mat, "mu_fixed_test")
    # pred_bias <- pred_res$pred_bias
    # pred_rmse <- pred_res$pred_rmse
    
  }
  
  # vars_to_round <- c("est_fwc", "se_fwc", "CI_l", "CI_u", "est_auc", "se_auc", "ci_auc", "pred_bias", "pred_rmse", "eps", "pred_full")
  vars_to_round <- c("est_fwc", "se_fwc", "CI_l", "CI_u", "est_auc", "se_auc", "ci_auc", "r0_sd", "r1_sd")
  for (i in seq_along(vars_to_round)) {
    assign(vars_to_round[i], round(get(vars_to_round[i]), 2))
  }
  summary_res <- c(
    sc00, 
    paste0(est_fwc, " (", se_fwc, ", CI: ", CI_l, ", ", CI_u, ")"), 
    paste0(est_auc, " (", se_auc, ", CI: ", ci_auc[1], ", ", ci_auc[2], ")"), 
    r0_sd, r1_sd, convergence
    # pred_bias, pred_rmse, eps, convergence, unlist(pred_full)
  ) 
  return(list(output = output, summary_res = summary_res))
}

### RUN Code ====================================
horseshoe_prior <- data.frame(
  lambda = c("student_t(1, 0, 1)", "student_t(1, 0, 1)", "student_t(1, 0, 1)", "student_t(1, 0, 1)", 
             rep("student_t(1, 0, 5)", 2)), 
  tau = c("student_t(1, 0, 1)", "student_t(1, 0, 2)", "student_t(1, 0, 5)", "student_t(1, 0, 10)", 
          "student_t(1, 0, 5)", "student_t(1, 0, 10)")
)
ig_prior <-  data.frame(
  a_slab = c(0.1, 0.001), b_slab = c(0.1, 0.001)
)
all.config <- bind_rows(
  # LMM
  expand.grid(
    trt.group = c(1, 0),
    model.type = c(1:3),
    lmm.type = c(1, 2), 
    gmm.type = 1,
    k.knots = c(4),
    data.type = c(1:2),
    sigma_prior = c("student_t(1, 0, 5)"), 
    hs_prior = 1,
    ig_prior = 1, 
    mcmc_seed = c(123)
  ), 
  # BLMM
  expand.grid(
    trt.group = c(1, 0),
    model.type = c(4:6),
    lmm.type = c(1, 2), 
    gmm.type = 1,
    k.knots = c(4),
    data.type = c(1:2),
    sigma_prior = c("student_t(1, 0, 5)"), 
    hs_prior = 1,
    ig_prior = 1, 
    mcmc_seed = c(123, 234, 345, 456, 567)
  ), 
  # sGMM
  expand.grid(
    trt.group = c(1, 0),
    model.type = 7,
    lmm.type = 1, 
    gmm.type = c(2, 4, 5),
    k.knots = c(4),
    data.type = c(1:2),
    sigma_prior = c("student_t(1, 0, 5)"),
    hs_prior = 1,
    ig_prior = 1, 
    mcmc_seed = c(123, 234, 345, 456, 567)
  ),
  # TL-sGMM
  # expand.grid(
  #   trt.group = c(1, 0),
  #   model.type = 7,
  #   lmm.type = 1, 
  #   gmm.type = c(6:10),
  #   k.knots = c(4),
  #   data.type = c(1:2),
  #   sigma_prior = c("student_t(1, 0, 5)"),
  #   hs_prior = c(1),
  #   ig_prior = 1, 
  #   mcmc_seed = c(123, 234, 345, 456, 567)
  # ),
  expand.grid(
    trt.group = c(1, 0),
    model.type = 7,
    lmm.type = 1, 
    gmm.type = c(11),
    k.knots = c(4),
    data.type = c(1:2),
    sigma_prior = c("student_t(1, 0, 5)"),
    hs_prior = c(1),
    ig_prior = c(2), 
    mcmc_seed = c(123, 234, 345, 456, 567)
  )
)

res <- main.func(
  data_all = data_all,
  trt.group = all.config$trt.group[sc00],
  n_cls = 5, 
  var_inflate = 5, 
  a0 = 10, 
  model.type = all.config$model.type[sc00],
  lmm.type = all.config$lmm.type[sc00],
  gmm.type = all.config$gmm.type[sc00], 
  k1 = all.config$k.knots[sc00], 
  k2 = 2, 
  data.type = all.config$data.type[sc00],
  sigma_prior = as.character(all.config$sigma_prior)[sc00], 
  lambda_prior = as.character(horseshoe_prior$lambda[all.config$hs_prior[sc00]]), 
  tau_prior = as.character(horseshoe_prior$tau[all.config$hs_prior[sc00]]),
  a_slab = as.character(ig_prior$a_slab[all.config$ig_prior[sc00]]), 
  b_slab = as.character(ig_prior$b_slab[all.config$ig_prior[sc00]]), 
  mcmc_seed = all.config$mcmc_seed[sc00]
)
saveRDS(res$output, file = "./results/fit_sc00.RData")

### write file in a parallel way
resultpath <- "./results/output.csv"
lockfile <- "./results/lockfile.lock"
lock <- lock(lockfile, timeout = Inf)
on.exit(unlock(lock), add = TRUE)
if (!file.exists(resultpath)) {
  write_excel_csv(as.data.frame(t(res$summary_res)), file = resultpath, append = T, col_names = TRUE)
} else {
  write_excel_csv(as.data.frame(t(res$summary_res)), file = resultpath, append = T, col_names = FALSE)
}

unlock(lock)

