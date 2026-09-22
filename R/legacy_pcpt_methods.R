#' SWOG Analysis Pipeline
#'
#' This script defines functions for analyzing the SWOG dataset (including synthetic versions),
#' focusing on methods for handling linkage bias in survival analysis with time-dependent covariates.
#'
#' ## Core Components:
#' - **NLAC approach:** Implements the non-linkage adjusted censoring method under time-dependent Cox models
#' - **IPW approach:** Uses inverse probability weighting to account for non-random linkage
#' - **Complete Case (CC) **
#' - **Augmented Complete Case (ACC) methods**
#'
#' ## Functions:
#' - `tdcSWOG_nlac_analysis_robust()`: NLAC Cox model with two change points (e.g., at 6.5 and 7.5 years)
#' - `tdc_real_data_analysis()`: Applies all four methods (IPW, CC, ACC, NLAC), summarizes estimates, SEs,
#'   and confidence intervals for treatment effect changes at both change points.
#' - `ti_real_data_analysis()`: Applies all four methods (IPW, CC, ACC, NLAC) using only time
#'   independent covariates, summarizes estimates, SEs,
#'   and confidence intervals for treatment effect.
#'
#' ## Dependencies:
#' - `survival`: Cox regression and survival object creation
#' - `dplyr`: Data manipulation
#' - `multcomp`: Inference for linear combinations (e.g., treatment effect after change points)
#'
#' ## Input Data Requirements:
#' - `data` is a list with:
#'   - `df`: per-subject data frame with columns `id`, `Y`, `Delta`, `linkage`, `txg`, `Y_1`, `Delta_1`, etc.
#'   - `tdata`: time-dependent data generated via `tmerge()`
#'
#'
expit <- function(x) 1 / (1 + exp(-x))
ci_95_constant = qnorm(0.975)


#' @param data A list with:
#'   - `df`: a data frame with one row per subject, containing:
#'     - `id`: subject ID
#'     - `Y`, `Delta`: observed time and event indicator (linked)
#'     - `Y_1`, `Delta_1`: observed time and indicator from the trial-only period
#'     - `linkage`: binary indicator for linkage status
#'     - `txg`: treatment indicator
#'   - `tdata`: not used directly (recreated from `df` internally)
#'
#' @param cox_covariates A character vector of covariate names (e.g., `c("txg", "T1", "T2")`) to include in the model.
#'
#' @return A `coxph` object with robust variance estimates (`cluster(id)`).
#'
tdcSWOG_nlac_analysis_robust <- function(data, cox_covariates) {
  df = data$df
  tdata <- data$tdata
  n = nrow(df)
  linked = which(df$linkage == 1)

  # Update Y and Delta for unlinked people
  df$Y[-linked] = df$Y_1[-linked]
  df$Delta[-linked] = df$Delta_1[-linked]

  # Now add time-dependent covariates.
  first_change_point = 2373 # first change point at 6.5 years.
  # first_change_point = 19 # for testing code
  t_1 = rep(first_change_point, n)
  tdata = tmerge(df, df, id=id, death = event(Y, Delta),
                 T1 = tdc(t_1, txg))
  tdata$T1[is.na(tdata$T1)] = 0

  second_change_point = 2738 # second change point at 7.5 years.
  # second_change_point = 21 # for testing code
  t_2 = rep(second_change_point, n)
  tdata = tmerge(tdata, df, id=id, T2 = tdc(t_2, txg))
  tdata$T2[is.na(tdata$T2)] = 0

  # --- Cox model design
  Z_cox_formula <- paste(cox_covariates, collapse = " + ")
  formula_str <- as.formula(paste("Surv(tstart, tstop, death) ~", Z_cox_formula,
                                  ' + cluster(id)'))
  res = coxph(formula_str, data=tdata)
  return(res)
}

#' Run All 4 methods (IPW, CC, ACC, NLAC)  on SWOG Data
#'
#' Returns coefficient estimates, standard errors, and confidence intervals for
#' treatment effect at post-change time points.
#'
#' @param data A list containing:
#'   - `df`: per-subject data frame as described in `tdcSWOG_nlac_analysis_robust`
#'   - `tdata`: counting process format data frame (optional)
#'
#' @param cv A character vector of covariate names for the Cox models (e.g., `c("txg", "T1", "T2")`)
#' @param lv A character vector of covariate names used for modeling linkage in the IPW analysis
#'
#' @return A list containing:
#'   - `beta_hat`: matrix of estimates (rows: methods, cols: covariates)
#'   - `beta_hat_se`: matrix of standard errors
#'   - `T1_conf`: 95% CI for treatment effect after first change point
#'   - `T2_conf`: 95% CI for treatment effect after second change point
#'   - `result_ipw`, `result_cc`, `result_acc`, `result_nlac`: full model objects from each method
#'
#' @details
#' This function computes:
#' - treatment effect estimates from all four methods
#' - linear combinations such as `txg + T1` and `txg + T1 + T2`
#' - confidence intervals using both analytic (IPW) and multcomp-based (others) approaches

tdc_real_data_analysis <- function(data, cv=td_cox_covariates,
                               lv=c(cox_covariates, additional_logistic_covariates)) {

  result_ipw =  tdc_ipw_analysis(data, cox_covariates = cv,
                                 lv)
  beta_hat_ipw = result_ipw$beta_hat
  beta_hat_ipw_se = result_ipw$beta_hat_se
  names(beta_hat_ipw_se) = names(beta_hat_ipw)

  result_cc = tdc_complete_case_analysis_robust(data,
                                                cox_covariates = cv)
  beta_hat_cc = result_cc$coefficients
  beta_hat_cc_se = sqrt(diag(result_cc$var))
  names(beta_hat_cc_se) = names(beta_hat_cc)

  result_acc = tdc_augmented_complete_case_analysis_robust(data,
                                                           cox_covariates = cv)
  beta_hat_acc = result_acc$coefficients
  beta_hat_acc_se = sqrt(diag(result_acc$var))
  names(beta_hat_acc_se) = names(beta_hat_acc)

  result_nlac = tdcSWOG_nlac_analysis_robust(data, cv)
  beta_hat_nlac = result_nlac$coefficients
  beta_hat_nlac_se = sqrt(diag(result_nlac$var))
  names(beta_hat_nlac_se) = names(beta_hat_nlac)

  beta_hat = rbind(beta_hat_ipw, beta_hat_cc, beta_hat_acc, beta_hat_nlac)
  beta_hat_se = rbind(beta_hat_ipw_se,
                      beta_hat_cc_se,
                      beta_hat_acc_se,
                      beta_hat_nlac_se)
  rownames(beta_hat) = c("IPW", "CC", "ACC", "NLAC")
  rownames(beta_hat_se) = c("IPW", "CC", "ACC", "NLAC")

  # Added some code to compute confidence intervals
  ## For ipw
  n = nrow(data$df)
  features = names(result_ipw$beta_hat)
  beta_var = result_ipw$beta_var
  txg_loc = which(features == "txg")
  T1_loc = which(features == "T1")
  after_T1_estimate = result_ipw$beta_hat["txg"] + result_ipw$beta_hat["T1"]

  after_T1_cov = beta_var[txg_loc, txg_loc] + beta_var[T1_loc, T1_loc] + 2 * beta_var[txg_loc, T1_loc]
  after_T1_se = sqrt(after_T1_cov / n)

  result_ipw$after_T1_cf_lower = after_T1_estimate - ci_95_constant * after_T1_se
  result_ipw$after_T1_cf_upper = after_T1_estimate + ci_95_constant * after_T1_se

  T2_loc = which(features == "T2")
  after_T2_estimate = result_ipw$beta_hat["txg"] + result_ipw$beta_hat["T1"] +
    result_ipw$beta_hat["T2"]
  after_T2_cov = beta_var[txg_loc, txg_loc] + beta_var[T1_loc, T1_loc] + beta_var[T2_loc, T2_loc] +
    2 * (beta_var[txg_loc, T1_loc] + beta_var[txg_loc, T2_loc] + beta_var[T1_loc, T2_loc])
  after_T2_se = sqrt(after_T2_cov / n)
  result_ipw$after_T2_cf_lower = after_T2_estimate - ci_95_constant * after_T2_se
  result_ipw$after_T2_cf_upper = after_T2_estimate + ci_95_constant * after_T2_se

  T1_conf = matrix(0, 4, 2)  # since we have four methods
  # compute the confidence intervals for parameter estimate after 7.5 years
  conf_int_T1_cc = confint(glht(result_cc, linfct="txg+T1==0"))
  conf_int_T1_acc = confint(glht(result_acc, linfct="txg+T1==0"))
  conf_int_T1_nlac = confint(glht(result_nlac, linfct="txg+T1==0"))
  T1_conf[1, 1] = result_ipw$after_T1_cf_lower
  T1_conf[1, 2] = result_ipw$after_T1_cf_upper
  T1_conf[2, 1] = conf_int_T1_cc$confint[2]
  T1_conf[2, 2] = conf_int_T1_cc$confint[3]
  T1_conf[3, 1] = conf_int_T1_acc$confint[2]
  T1_conf[3, 2] = conf_int_T1_acc$confint[3]
  T1_conf[4, 1] = conf_int_T1_nlac$confint[2]
  T1_conf[4, 2] = conf_int_T1_nlac$confint[3]
  rownames(T1_conf) = c("IPW", "CC", "ACC", "NLAC")

  T2_conf = matrix(0, 4, 2)

  conf_int_T2_cc = confint(glht(result_cc, linfct="txg+T1+T2==0"))
  conf_int_T2_acc = confint(glht(result_acc, linfct="txg+T1+T2==0"))
  conf_int_T2_nlac = confint(glht(result_nlac, linfct="txg+T1+T2==0"))
  T2_conf[1, 1] = result_ipw$after_T2_cf_lower
  T2_conf[1, 2] = result_ipw$after_T2_cf_upper
  T2_conf[2, 1] = conf_int_T2_cc$confint[2]
  T2_conf[2, 2] = conf_int_T2_cc$confint[3]
  T2_conf[3, 1] = conf_int_T2_acc$confint[2]
  T2_conf[3, 2] = conf_int_T2_acc$confint[3]
  T2_conf[4, 1] = conf_int_T2_nlac$confint[2]
  T2_conf[4, 2] = conf_int_T2_nlac$confint[3]
  rownames(T2_conf) = c("IPW", "CC", "ACC", "NLAC")

  return(list("beta_hat" = beta_hat,
              "beta_hat_se" = beta_hat_se,
              "T1_conf" = T1_conf,
              "T2_conf" = T2_conf,
              "result_ipw"=result_ipw,
              "result_cc"=result_cc,
              "result_acc"=result_acc,
              "result_nlac"=result_nlac))
}


#' Run All 4 methods (IPW, CC, ACC, NLAC)  on SWOG Data
#'
#' Returns coefficient estimates, standard errors, and confidence intervals for
#' treatment effect at post-change time points.
#'
#' @param data A list containing:
#'   - `df`: per-subject data frame as described in `tdcSWOG_nlac_analysis_robust`
#'   - `tdata`: counting process format data frame (optional)
#'
#' @param cv A character vector of covariate names for the Cox models (e.g., `c("txg", "T1", "T2")`)
#' @param lv A character vector of covariate names used for modeling linkage in the IPW analysis
#'
#' @return A list containing:
#'   - `beta_hat`: matrix of estimates (rows: methods, cols: covariates)
#'   - `beta_hat_se`: matrix of standard errors
#'   - `T1_conf`: 95% CI for treatment effect after first change point
#'   - `T2_conf`: 95% CI for treatment effect after second change point
#'   - `result_ipw`, `result_cc`, `result_acc`, `result_nlac`: full model objects from each method
#'
#' @details
#' This function computes:
#' - treatment effect estimates from all four methods
#' - linear combinations such as `txg + T1` and `txg + T1 + T2`
#' - confidence intervals using both analytic (IPW) and multcomp-based (others) approaches

ti_real_data_analysis <- function(data, cv=td_cox_covariates,
                                   lv=c(cox_covariates, additional_logistic_covariates)) {

  result_ipw =  ti_ipw_analysis(data, cox_covariates = cv,
                                 lv)
  beta_hat_ipw = result_ipw$beta_hat
  beta_hat_ipw_se = result_ipw$beta_hat_se
  names(beta_hat_ipw_se) = names(beta_hat_ipw)

  result_cc = ti_complete_case_analysis(data,
                                                cox_covariates = cv)
  beta_hat_cc = result_cc$coefficients
  beta_hat_cc_se = sqrt(diag(result_cc$var))
  names(beta_hat_cc_se) = names(beta_hat_cc)

  result_acc = ti_augmented_complete_case_analysis(data,
                                                           cox_covariates = cv)
  beta_hat_acc = result_acc$coefficients
  beta_hat_acc_se = sqrt(diag(result_acc$var))
  names(beta_hat_acc_se) = names(beta_hat_acc)

  result_nlac = ti_nlac_analysis(data, cv)
  beta_hat_nlac = result_nlac$coefficients
  beta_hat_nlac_se = sqrt(diag(result_nlac$var))
  names(beta_hat_nlac_se) = names(beta_hat_nlac)

  beta_hat = rbind(beta_hat_ipw, beta_hat_cc, beta_hat_acc, beta_hat_nlac)
  beta_hat_se = rbind(beta_hat_ipw_se,
                      beta_hat_cc_se,
                      beta_hat_acc_se,
                      beta_hat_nlac_se)
  rownames(beta_hat) = c("IPW", "CC", "ACC", "NLAC")
  rownames(beta_hat_se) = c("IPW", "CC", "ACC", "NLAC")

  return(list("beta_hat" = beta_hat,
              "beta_hat_se" = beta_hat_se,
              "result_ipw"=result_ipw,
              "result_cc"=result_cc,
              "result_acc"=result_acc,
              "result_nlac"=result_nlac))
}
