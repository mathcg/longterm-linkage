# linkage_functions.R
# ----------------------------------
# This script defines all functions used for simulation and estimation in 
# paper "Long-term effect estimation when combining clinical trial and observational follow-up
# datasets" by Cheng et al
#
# Contents:
# - Data generation functions for simulation settings in the manuscript
# - Estimation functions (Full cohort, Complete case, Augmented complete case, NLAC, IPW, AIPW, etc.)
# - Utility functions for risk set, gradient, and variance estimation
# - Wrapper function: `all_simulation()` to simulate and analyze one run
#
# Usage:
# This file is intended to be sourced from a wrapper script (e.g., run_all.R).
# It does not execute any simulations directly on load.
#
# Author: Gang Cheng
# Updated: April 2, 2025 by Yingqi Zhao
# ----------------------------------

expit <- function(x) 1 / (1 + exp(-x))

library(survival)
library(MASS)
library(dplyr)
library(collapse)

## time dependent setting, main paper
tdc_data_generation_setting_1 <- function(n, seed) {
  set.seed(seed)
  
  # True parameters
  beta_1 <- -log(4)
  beta_2 <- log(1.5)
  beta_3 <- 0.5
  lambda <- 0.06
  tau <- 5
  tau_M <- 16
  
  # Covariates
  Z_1 <- rbinom(n, 1, 0.5)             # binary treatment
  Z_2 <- rnorm(n, mean = 1, sd = 1)    # continuous covariate
  
  # Simulate survival time using piecewise exponential model
  t_0 <- rep(tau, n)
  U <- runif(n)
  log_U <- -log(1 - U)
  
  temp_pre <- lambda * exp(beta_1 * Z_1 + beta_2 * Z_2) * t_0
  T_failure <- numeric(n)
  fail_before_tau <- log_U < temp_pre
  
  # Case 1: fail before tau
  T_failure[fail_before_tau] <- log_U[fail_before_tau] * t_0[fail_before_tau] / temp_pre[fail_before_tau]
  
  # Case 2: fail after tau
  temp_post <- lambda * exp(beta_1 * Z_1[!fail_before_tau] +
                              beta_2 * Z_2[!fail_before_tau] +
                              beta_3 * Z_1[!fail_before_tau]) * t_0[!fail_before_tau]
  T_failure[!fail_before_tau] <- (log_U[!fail_before_tau] - temp_pre[!fail_before_tau] + temp_post) / 
    (temp_post / t_0[!fail_before_tau])
  
  #  Generate censoring times
  C_1 <- pmin(rexp(n, rate = 0.01 * Z_1 + 0.03), tau)
  tau_start <- C_1
  C_2 <- pmin(rexp(n, rate = 0.05 * Z_1 + 0.03) + tau_start, tau_M)
  C_time <- pmax(C_1, C_2) 
  
  #  Compute observed time and event indicator
  Y <- pmin(T_failure, C_time)
  Delta <- as.numeric(T_failure <= C_time)
  censored_in_trial <- as.numeric(C_1 < T_failure)
  failed_in_trial <- as.numeric(T_failure <= C_1)
  
  Y_1 <- pmin(T_failure, C_1) 
  
  # Generating linkage data
  linkage = rbinom(n, 1, prob=0.5)
  
  # Build initial dataset
  data <- data.frame(id = 1:n, Z_1 = Z_1, Z_2 = Z_2, Y = Y, Delta = Delta,
                     linkage = linkage, censored_in_trial = censored_in_trial,
                     failed_in_trial = failed_in_trial, Y_1 = Y_1,
                     "t_0"=t_0)
  
  # Create time-dependent covariate Z_3 = Z_1 * I(t ≥ τ)
  tdata <- tmerge(data, data, id = id, death = event(Y, Delta),
                  Z_3=tdc(t_0, Z_1))
  tdata$Z_3[is.na(tdata$Z_3)] <- 0 
   
  return(list("tdata"=tdata, "df"=data))
}



## Setting 2, CLAR case

tdc_data_generation_setting_2 <- function(n, seed) {
  set.seed(seed)
  
  # True parameters
  beta_1 <- -log(4)
  beta_2 <- log(1.5)
  beta_3 <- 0.5
  lambda <- 0.06
  tau <- 5
  tau_M <- 16
  
  # Covariates
  Z_1 <- rbinom(n, 1, 0.5)             # binary treatment
  Z_2 <- rnorm(n, mean = 1, sd = 1)    # continuous covariate
  
  # Simulate survival time using piecewise exponential model
  t_0 <- rep(tau, n)
  U <- runif(n)
  log_U <- -log(1 - U)
  
  temp_pre <- lambda * exp(beta_1 * Z_1 + beta_2 * Z_2) * t_0
  T_failure <- numeric(n)
  fail_before_tau <- log_U < temp_pre
  
  # Case 1: fail before tau
  T_failure[fail_before_tau] <- log_U[fail_before_tau] * t_0[fail_before_tau] / temp_pre[fail_before_tau]
  
  # Case 2: fail after tau
  temp_post <- lambda * exp(beta_1 * Z_1[!fail_before_tau] +
                              beta_2 * Z_2[!fail_before_tau] +
                              beta_3 * Z_1[!fail_before_tau]) * t_0[!fail_before_tau]
  T_failure[!fail_before_tau] <- (log_U[!fail_before_tau] - temp_pre[!fail_before_tau] + temp_post) / 
    (temp_post / t_0[!fail_before_tau])
  
  #  Generate censoring times
  C_1 <- pmin(rexp(n, rate = 0.01 * Z_1 + 0.03), tau)
  tau_start <- C_1
  C_2 <- pmin(rexp(n, rate = 0.05 * Z_1 + 0.03) + tau_start, tau_M)
  C_time <- pmax(C_1, C_2) 
  
  #  Compute observed time and event indicator
  Y <- pmin(T_failure, C_time)
  Delta <- as.numeric(T_failure <= C_time)
  censored_in_trial <- as.numeric(C_1 < T_failure)
  failed_in_trial <- as.numeric(T_failure <= C_1)
  
  Y_1 <- pmin(T_failure, C_1) 
  
  # Generating linkage data
  censored_in_trial_idx = which(C_1 < T_failure)
  gamma_1 = 0.5
  gamma_2 = 0.5
  
  linkage = rep(0, n)
  linkage[censored_in_trial_idx] = rbinom(length(censored_in_trial_idx), 1, 
                                      prob=expit(-0.25 + gamma_1 * Z_1[censored_in_trial_idx] + 
                                                   gamma_2 * Z_2[censored_in_trial_idx]))
  linkage[-censored_in_trial_idx] = rbinom(n - length(censored_in_trial_idx), 1, prob = 0.5)
  
  
  # Build initial dataset
  data <- data.frame(id = 1:n, Z_1 = Z_1, Z_2 = Z_2, Y = Y, Delta = Delta,
                     linkage = linkage, censored_in_trial = censored_in_trial,
                     failed_in_trial = failed_in_trial, Y_1 = Y_1,
                     "t_0"=t_0)
  
  # Create time-dependent covariate Z_3 = Z_1 * I(t ≥ τ)
  tdata <- tmerge(data, data, id = id, death = event(Y, Delta),
                  Z_3=tdc(t_0, Z_1))
  tdata$Z_3[is.na(tdata$Z_3)] <- 0 
  
  return(list("tdata"=tdata, "df"=data))
}


## Setting 3, LNAR(\tilde T) case
tdc_data_generation_setting_3 <- function(n, seed) {
  set.seed(seed)
  
  # True parameters
  beta_1 <- -log(4)
  beta_2 <- log(1.5)
  beta_3 <- 0.5
  lambda <- 0.06
  tau <- 5
  tau_M <- 16
  
  # Covariates
  Z_1 <- rbinom(n, 1, 0.5)             # binary treatment
  Z_2 <- rnorm(n, mean = 1, sd = 1)    # continuous covariate
  
  # Simulate survival time using piecewise exponential model
  t_0 <- rep(tau, n)
  U <- runif(n)
  log_U <- -log(1 - U)
  
  temp_pre <- lambda * exp(beta_1 * Z_1 + beta_2 * Z_2) * t_0
  T_failure <- numeric(n)
  fail_before_tau <- log_U < temp_pre
  
  # Case 1: fail before tau
  T_failure[fail_before_tau] <- log_U[fail_before_tau] * t_0[fail_before_tau] / temp_pre[fail_before_tau]
  
  # Case 2: fail after tau
  temp_post <- lambda * exp(beta_1 * Z_1[!fail_before_tau] +
                              beta_2 * Z_2[!fail_before_tau] +
                              beta_3 * Z_1[!fail_before_tau]) * t_0[!fail_before_tau]
  T_failure[!fail_before_tau] <- (log_U[!fail_before_tau] - temp_pre[!fail_before_tau] + temp_post) / 
    (temp_post / t_0[!fail_before_tau])
  
  #  Generate censoring times
  C_1 <- pmin(rexp(n, rate = 0.01 * Z_1 + 0.03), tau)
  tau_start <- C_1
  C_2 <- pmin(rexp(n, rate = 0.05 * Z_1 + 0.03) + tau_start, tau_M)
  C_time <- pmax(C_1, C_2) 
  
  #  Compute observed time and event indicator
  Y <- pmin(T_failure, C_time)
  Delta <- as.numeric(T_failure <= C_time)
  censored_in_trial <- as.numeric(C_1 < T_failure)
  failed_in_trial <- as.numeric(T_failure <= C_1) 
  
  Y_1 <- pmin(T_failure, C_1) 

  
  # Generating linkage data
  censored_in_trial_idx = which(C_1 < T_failure)
  gamma_1 = 0.5
  gamma_2 = 0.5
  
  linkage = rep(0, n)
  linkage[censored_in_trial_idx] = rbinom(length(censored_in_trial_idx), 1, 
                                      prob=expit(-0.25 + gamma_1 * Z_1[censored_in_trial_idx] + 
                                                   gamma_2 * Z_2[censored_in_trial_idx]
                                                 -0.01 * T_failure[censored_in_trial_idx] - 
                                                   0.01 * Delta[censored_in_trial_idx]))
  linkage[-censored_in_trial_idx] = rbinom(n - length(censored_in_trial_idx), 1, prob = 0.5)
  
  # Build initial dataset
  data <- data.frame(id = 1:n, Z_1 = Z_1, Z_2 = Z_2, Y = Y, Delta = Delta,
                     linkage = linkage, censored_in_trial = censored_in_trial,
                     failed_in_trial = failed_in_trial, Y_1 = Y_1,
                     "t_0"=t_0)
  
  
  # Create time-dependent covariate Z_3  = Z_1 * I(t ≥ τ)
  tdata <- tmerge(data, data, id = id, death = event(Y, Delta),
                  Z_3=tdc(t_0, Z_1))
  tdata$Z_3[is.na(tdata$Z_3)] <- 0 
  
  return(list("tdata"=tdata, "df"=data))
}


## Setting 4, LNAR(C_2) case
tdc_data_generation_setting_4 <- function(n, seed) {
  set.seed(seed)
  
  # True parameters
  beta_1 <- -log(4)
  beta_2 <- log(1.5)
  beta_3 <- 0.5
  lambda <- 0.06
  tau <- 5
  tau_M <- 16
  
  # Covariates
  Z_1 <- rbinom(n, 1, 0.5)             # binary treatment
  Z_2 <- rnorm(n, mean = 1, sd = 1)    # continuous covariate
  
  # Simulate survival time using piecewise exponential model
  t_0 <- rep(tau, n)
  U <- runif(n)
  log_U <- -log(1 - U)
  
  temp_pre <- lambda * exp(beta_1 * Z_1 + beta_2 * Z_2) * t_0
  T_failure <- numeric(n)
  fail_before_tau <- log_U < temp_pre
  
  # Case 1: fail before tau
  T_failure[fail_before_tau] <- log_U[fail_before_tau] * t_0[fail_before_tau] / temp_pre[fail_before_tau]
  
  # Case 2: fail after tau
  temp_post <- lambda * exp(beta_1 * Z_1[!fail_before_tau] +
                              beta_2 * Z_2[!fail_before_tau] +
                              beta_3 * Z_1[!fail_before_tau]) * t_0[!fail_before_tau]
  T_failure[!fail_before_tau] <- (log_U[!fail_before_tau] - temp_pre[!fail_before_tau] + temp_post) / 
    (temp_post / t_0[!fail_before_tau])
  
  #  Generate censoring times
  C_1 <- pmin(rexp(n, rate = 0.01 * Z_1 + 0.03), tau)
  tau_start <- C_1
  C_2 <- pmin(rexp(n, rate = 0.05 * Z_1 + 0.03) + tau_start, tau_M)
  C_time <- pmax(C_1, C_2) 
  
  #  Compute observed time and event indicator
  Y <- pmin(T_failure, C_time)
  Delta <- as.numeric(T_failure <= C_time)
  censored_in_trial <- as.numeric(C_1 < T_failure)
  failed_in_trial <- as.numeric(T_failure <= C_1)  
  
  Y_1 <- pmin(T_failure, C_1) 
  
  
  # Generating linkage data
  censored_in_trial_idx = which(C_1 < T_failure)
  gamma_1 = 0.5
  gamma_2 = 0.5
  
  linkage = rep(0, n)
  linkage[censored_in_trial_idx] = rbinom(length(censored_in_trial_idx), 1, 
                                      prob=expit(-0.25 + gamma_1 * Z_1[censored_in_trial_idx] + 
                                                   gamma_2 * Z_2[censored_in_trial_idx]
                                                 -0.1 * C_2[censored_in_trial_idx] - 
                                                   0.1 * Delta[censored_in_trial_idx]))
  linkage[-censored_in_trial_idx] = rbinom(n - length(censored_in_trial_idx), 1, prob = 0.5)
  
  # Build initial dataset
  data <- data.frame(id = 1:n, Z_1 = Z_1, Z_2 = Z_2, Y = Y, Delta = Delta,
                     linkage = linkage, censored_in_trial = censored_in_trial,
                     failed_in_trial = failed_in_trial,  Y_1 = Y_1,
                     "t_0"=t_0)
  
  # Create time-dependent covariate Z_3  = Z_1 * I(t ≥ τ)
  tdata <- tmerge(data, data, id = id, death = event(Y, Delta),
                  Z_3=tdc(t_0, Z_1))
  tdata$Z_3[is.na(tdata$Z_3)] <- 0 
  
  return(list("tdata"=tdata, "df"=data))
}


## Setting 5, MAR case
tdc_data_generation_setting_5 <- function(n, seed) {
  set.seed(seed)
  
  # True parameters
  beta_1 <- -log(4)
  beta_2 <- log(1.5)
  beta_3 <- 0.5
  lambda <- 0.06
  tau <- 5
  tau_M <- 16
  
  # Covariates
  Z_1 <- rbinom(n, 1, 0.5)             # binary treatment
  Z_2 <- rnorm(n, mean = 1, sd = 1)    # continuous covariate
  
  # Simulate survival time using piecewise exponential model
  t_0 <- rep(tau, n)
  U <- runif(n)
  log_U <- -log(1 - U)
  
  temp_pre <- lambda * exp(beta_1 * Z_1 + beta_2 * Z_2) * t_0
  T_failure <- numeric(n)
  fail_before_tau <- log_U < temp_pre
  
  # Case 1: fail before tau
  T_failure[fail_before_tau] <- log_U[fail_before_tau] * t_0[fail_before_tau] / temp_pre[fail_before_tau]
  
  # Case 2: fail after tau
  temp_post <- lambda * exp(beta_1 * Z_1[!fail_before_tau] +
                              beta_2 * Z_2[!fail_before_tau] +
                              beta_3 * Z_1[!fail_before_tau]) * t_0[!fail_before_tau]
  T_failure[!fail_before_tau] <- (log_U[!fail_before_tau] - temp_pre[!fail_before_tau] + temp_post) / 
    (temp_post / t_0[!fail_before_tau])
  
  #  Generate censoring times
  C_1 <- pmin(rexp(n, rate = 0.01 * Z_1 + 0.03), tau)
  tau_start <- C_1
  C_2 <- pmin(rexp(n, rate = 0.05 * Z_1 + 0.03) + tau_start, tau_M)
  C_time <- pmax(C_1, C_2) 
  
  #  Compute observed time and event indicator
  Y <- pmin(T_failure, C_time)
  Delta <- as.numeric(T_failure <= C_time)
  censored_in_trial <- as.numeric(C_1 < T_failure)
  failed_in_trial <- as.numeric(T_failure <= C_1)  
  
  Y_1 <- pmin(T_failure, C_1) 
  
  # Generating linkage data 
  gamma_1 = 0.5
  gamma_2 = 0.5
  
  linkage = rep(0, n)
  linkage = rbinom(n, 1, prob=expit(-0.25 + gamma_1 * Z_1 + gamma_2 * Z_2)) 
  
  
  # Build initial dataset
  data <- data.frame(id = 1:n, Z_1 = Z_1, Z_2 = Z_2, Y = Y, Delta = Delta,
                     linkage = linkage, censored_in_trial = censored_in_trial,
                     failed_in_trial = failed_in_trial,  Y_1 = Y_1,
                     "t_0"=t_0)
  
  # Create time-dependent covariate Z_3  = Z_1 * I(t ≥ τ)
  tdata <- tmerge(data, data, id = id, death = event(Y, Delta),
                  Z_3=tdc(t_0, Z_1))
  tdata$Z_3[is.na(tdata$Z_3)] <- 0 
  
  return(list("tdata"=tdata, "df"=data))
}


## supplementary H.3, linkage assumption satisfies MAR  
ti_data_generation_setting_6 <- function(n, seed) {
  set.seed(seed)
  
  # True parameters
  beta_1 = -log(2)
  beta_2 = log(2)
  beta_3 = 0.2
  lambda <- 0.06
  tau <- 7
  tau_M <- 16
  
  # Covariates
  Z_1 <- rbinom(n, 1, 0.5)             # binary treatment
  Z_2 = rnorm(n, mean=-1, sd=1)
  Z_3 = rnorm(n, mean=1, sd=2)        # continuous covariate
  
  # Simulate survival time using piecewise exponential model
  beta_true = c(beta_1, beta_2, beta_3)
  
  lambda = 0.15 # assume that the hazard function is constant, correspond to exponential distribution.
  T_uniform = runif(n)
  T_failure = -log(1 - T_uniform) / (lambda * exp(beta_1 * Z_1 + beta_2 * Z_2 + beta_3 * Z_3)) # lambda is a constant
  
  
  #  Generate censoring times
  C_1 <- runif(n, min = 0, max = 3.5)
  tau_start <- C_1
  C_2 <- pmin(rexp(n, rate = 0.8 * Z_1 + 0.03) + 1, tau_M)
  C_time <- pmax(C_1, C_2) 
  
  #  Compute observed time and event indicator
  Y <- pmin(T_failure, C_time)
  Delta <- as.numeric(T_failure <= C_time)
  censored_in_trial <- as.numeric(C_1 < T_failure)
  failed_in_trial <- as.numeric(T_failure <= C_1)  
  
  Y_1 <- pmin(T_failure, C_1) 
  
  # Generating linkage data 
  gamma_1 = 0.5
  gamma_2 = 0.5
  
  linkage = rep(0, n)
  linkage = rbinom(n, 1, prob=expit(-0.25 + gamma_1 * Z_1 + gamma_2 * Z_2)) 
  
  
  # Build initial dataset
  data <- data.frame(id = 1:n, Z_1 = Z_1, Z_2 = Z_2, Z_3 = Z_3, 
                     Y = Y, Delta = Delta, linkage = linkage,
                     censored_in_trial = censored_in_trial,
                     failed_in_trial = failed_in_trial,  Y_1 = Y_1)
  
  return(list("df"= data))
}




## time independent setting, supplementary H.3, linkage satisfies CLAR
ti_data_generation_setting_7 <- function(n, seed) {
  set.seed(seed)
  
  # True parameters
  beta_1 = -log(2)
  beta_2 = log(2)
  beta_3 = 0.2
  lambda <- 0.06
  tau <- 7
  tau_M <- 16
  
  # Covariates
  Z_1 <- rbinom(n, 1, 0.5)             # binary treatment
  Z_2 = rnorm(n, mean=-1, sd=1)
  Z_3 = rnorm(n, mean=1, sd=2)        # continuous covariate
  
  # Simulate survival time using piecewise exponential model
  beta_true = c(beta_1, beta_2, beta_3)
  
  lambda = 0.15 # assume that the hazard function is constant, correspond to exponential distribution.
  T_uniform = runif(n)
  T_failure = -log(1 - T_uniform) / (lambda * exp(beta_1 * Z_1 + beta_2 * Z_2 + beta_3 * Z_3)) # lambda is a constant
  
  
  #  Generate censoring times
  C_1 <- runif(n, min = 0, max = 3.5)
  tau_start <- C_1
  C_2 <- pmin(rexp(n, rate = 0.8 * Z_1 + 0.03) + 1, tau_M)
  C_time <- pmax(C_1, C_2) 
  
  #  Compute observed time and event indicator
  Y <- pmin(T_failure, C_time)
  Delta <- as.numeric(T_failure <= C_time)
  censored_in_trial <- as.numeric(C_1 < T_failure)
  failed_in_trial <- as.numeric(T_failure <= C_1) 
  
  Y_1 <- pmin(T_failure, C_1) 
  
  # Generating linkage data
  # Assumed that P(L = 1 | C_1 < T_failure, Z) = expit(gamma^T Z)
  # feels to me that in simulation, I should simply generate data
  # such that P(L = 1 | Z) = expit(gamma^T Z). 
  # after_tau = which(C_1 >= tau & T_failure >= tau)
  
  censored_in_trial_idx = which(C_1 < T_failure)
  gamma_1 = 0.5
  gamma_2 = 0.5
  
  linkage = rep(0, n)
  linkage[censored_in_trial_idx] = rbinom(length(censored_in_trial_idx), 1, 
                                      prob=expit(-0.25 + gamma_1 * Z_1[censored_in_trial_idx] + 
                                                   gamma_2 * Z_2[censored_in_trial_idx]))
  linkage[-censored_in_trial_idx] = rbinom(n - length(censored_in_trial_idx), 1, prob = 0.5)
   
  # Build initial dataset
  data <- data.frame(id = 1:n, Z_1 = Z_1, Z_2 = Z_2, Z_3 = Z_3,
                     Y = Y, Delta = Delta,
                     linkage = linkage, censored_in_trial = censored_in_trial,
                     failed_in_trial = failed_in_trial,  Y_1 = Y_1)
  
  return(list("df"= data))
}

## Below are analysis methods
## data is a list, containing a data frame from observed data, including  
## individual covariates such as Z_1, Z_2 etc
## linkage: linkage indicator
## Y: observed time 
## Delta: censoring indicator
## If data is time dependent, also includes a dataset with time-dependent 
## covariates, built using the tmerge() function from the survival package.


##############################################
##           Full cohort analysis           ##
##              (Oracle method)             ##
##                                          ##
##############################################


#' Time Independent Full Cohort Cox Analysis  
#'@param data A list containing:

#'   - `df`: a data frame with one row per subject and the following columns:
#'
#'     - `id`: unique subject identifier
#'     - `Z_1`, `Z_2`, ...: individual-level covariates
#'     - `Y`: observed time (either failure or censoring) assuming full linkage
#'     - `Delta`: failure indicator (1 = event, 0 = censored) in the linked dataset
#'     - `linkage`: linkage indicator (1 = linked to followup observational study, 0 = not linked)
#'     - `censored_in_trial`: binary indicator (1 = censored during trial period, 0 = otherwise)
#'     - `failed_in_trial`: binary indicator (1 = failed during trial period, 0 = otherwise)
#'     - `Y_1`: observed time in the trial-only data, assuming no linkage
#'     
#' @param cox_covariates A character vector of column names in `df` to be included as covariates in the Cox model.
#' @return A `coxph` model object fit to the linked full cohort data.

ti_full_cohort_analysis <- function(data, cox_covariates) {
  # For full cohort analysis,  everybody is linked to the medicare;
  # Further, we are making the "no gap" assumption. So, we can just compare
  # the survival time to the medicare censoring time.
  df = data$df
  n = nrow(df)
  
  # --- Cox model design
  Z_cox_formula <- paste(cox_covariates, collapse = " + ")
  formula_str <- as.formula(paste("Surv(Y, Delta) ~", Z_cox_formula)) 
  res = coxph(formula_str, data = df)
  
  return(res)
}

#' Time-Dependent Full Cohort Cox Analysis (Robust) 
#'
#' @param data A list containing:
#' #'   - `df`: a data frame with one row per subject and the following columns:
#'
#'     - `id`: unique subject identifier
#'     - `Z_1`, `Z_2`, ...: individual-level covariates
#'     - `Y`: observed time (either failure or censoring) assuming full linkage
#'     - `Delta`: failure indicator (1 = event, 0 = censored) in the linked dataset
#'     - `linkage`: linkage indicator (1 = linked to followup observational study, 0 = not linked)
#'     - `censored_in_trial`: binary indicator (1 = censored during trial period, 0 = otherwise)
#'     - `failed_in_trial`: binary indicator (1 = failed during trial period, 0 = otherwise)
#'     - `Y_1`: observed time in the trial-only data, assuming no linkage
#'     
#'   - `tdata`: a data frame in counting process format, with one or more rows per subject.
#'     It must include:
#'       - `id`: unique subject identifier
#'       - `tstart`, `tstop`: interval start and stop times
#'       - `death`: event indicator (1 = event, 0 = censored)
#'       - Covariates including any time-dependent ones such as `Z_3` 
#'       - t_0 (optional: change point in time
#'
#' @param cox_covariates A character vector of column names in `tdata` to be included as covariates in the Cox model.
#'
#' @return A `coxph` model object fit to the full cohort with robust variance via clustering by subject ID.
#' 
#'
tdc_full_cohort_analysis_robust <- function(data, cox_covariates) {
  tdata = data$tdata
  
  # --- Cox model design
  Z_cox_formula <- paste(cox_covariates, collapse = " + ")
  formula_str <- as.formula(paste("Surv(tstart, tstop, death) ~", Z_cox_formula,
                                  ' + cluster(id)'))
  res = coxph(formula_str, data = tdata) 
  return(res)
}



##############################################
##         Complete case analysis           ##
##                                          ##
##############################################


#' Time-Independent Complete Case Analysis 
#'
#'@param data A list containing:

#'   - `df`: a data frame with one row per subject and the following columns:
#'
#'     - `id`: unique subject identifier
#'     - `Z_1`, `Z_2`, ...: individual-level covariates
#'     - `Y`: observed time (either failure or censoring) assuming full linkage
#'     - `Delta`: failure indicator (1 = event, 0 = censored) in the linked dataset
#'     - `linkage`: linkage indicator (1 = linked to followup observational study, 0 = not linked)
#'     - `censored_in_trial`: binary indicator (1 = censored during trial period, 0 = otherwise)
#'     - `failed_in_trial`: binary indicator (1 = failed during trial period, 0 = otherwise)
#'     - `Y_1`: observed time in the trial-only data, assuming no linkage
#'     
#' @param cox_covariates A character vector of column names in `df` to be included as covariates in the Cox model.
#' @return A `coxph` model object fit to the linked (complete-case) data.
#'
ti_complete_case_analysis <- function(data, cox_covariates) {
  # Complete case only use patients that are linked.
  df = data$df
  n = nrow(df)
  linked_people = which(df$linkage == 1)
  data_cc = df[df$id %in% linked_people, ]
  
  # --- Cox model design
  Z_cox_formula <- paste(cox_covariates, collapse = " + ")
  formula_str <- as.formula(paste("Surv(Y, Delta) ~", Z_cox_formula)) 
  res = coxph(formula_str, data = data_cc) 
  
  return(res)
}



#' Time-Dependent Complete Case Analysis (Robust) 
#'
#' @param data A list containing:
#' #'   - `df`: a data frame with one row per subject and the following columns:
#'
#'     - `id`: unique subject identifier
#'     - `Z_1`, `Z_2`, ...: individual-level covariates
#'     - `Y`: observed time (either failure or censoring) assuming full linkage
#'     - `Delta`: failure indicator (1 = event, 0 = censored) in the linked dataset
#'     - `linkage`: linkage indicator (1 = linked to followup observational study, 0 = not linked)
#'     - `censored_in_trial`: binary indicator (1 = censored during trial period, 0 = otherwise)
#'     - `failed_in_trial`: binary indicator (1 = failed during trial period, 0 = otherwise)
#'     - `Y_1`: observed time in the trial-only data, assuming no linkage
#'     
#'   - `tdata`: a data frame in counting process format, with one or more rows per subject.
#'     It must include:
#'       - `id`: unique subject identifier
#'       - `tstart`, `tstop`: interval start and stop times
#'       - `death`: event indicator (1 = event, 0 = censored)
#'       - Covariates including any time-dependent ones such as `Z_3` 
#'       - t_0 (optional): change point in time
#'
#' @param cox_covariates A character vector of column names in `tdata` to be included as covariates in the Cox model.
#'
#' @return A `coxph` model object fit to the linked (complete-case) data with robust variance via clustering by subject ID.
#' 
tdc_complete_case_analysis_robust <- function(data, cox_covariates) {
  # Complete case only use patients that are linked.
  df = data$df
  n = nrow(df)
  linked_people = which(df$linkage == 1)
  tdata = data$tdata
  tdata_cc = tdata[tdata$id %in% linked_people, ]
  
  # --- Cox model design
  Z_cox_formula <- paste(cox_covariates, collapse = " + ")
  formula_str <- as.formula(paste("Surv(tstart, tstop, death) ~", Z_cox_formula,
                                  ' + cluster(id)'))
  res = coxph(formula_str, data = tdata_cc) 
  return(res)
}


##############################################
##                                          ##
##    Augmented Complete case analysis      ##
##                                          ##
##############################################

#' Time-Independent Augmented Complete Case Analysis  
#' 
#'@param data A list containing:

#'   - `df`: a data frame with one row per subject and the following columns:
#'
#'     - `id`: unique subject identifier
#'     - `Z_1`, `Z_2`, ...: individual-level covariates
#'     - `Y`: observed time (either failure or censoring) assuming full linkage
#'     - `Delta`: failure indicator (1 = event, 0 = censored) in the linked dataset
#'     - `linkage`: linkage indicator (1 = linked to followup observational study, 0 = not linked)
#'     - `censored_in_trial`: binary indicator (1 = censored during trial period, 0 = otherwise)
#'     - `failed_in_trial`: binary indicator (1 = failed during trial period, 0 = otherwise)
#'     - `Y_1`: observed time in the trial-only data, assuming no linkage
#'     
#' @param cox_covariates A character vector of column names in `df` to be included as covariates in the Cox model.
#' @return A `coxph` model object fit to the linked data, using the augmented complete case analysis.
#' 
## time independent
ti_augmented_complete_case_analysis <- function(data, cox_covariates) {
  # augmented complete case use patients that are linked or failed in trial.
  df = data$df
  n = nrow(df)
  linked= which(df$linkage == 1)
  failed_in_trial =  df$failed_in_trial  
  failed_in_trial_idx <- which(failed_in_trial == 1)
  linked_plus_failed_in_trial <- sort(unique(c(linked, failed_in_trial_idx))) 
  data_aug_cc = df[df$id %in% linked_plus_failed_in_trial , ]
  
  
  # --- Cox model design
  Z_cox_formula <- paste(cox_covariates, collapse = " + ")
  formula_str <- as.formula(paste("Surv(Y, Delta) ~", Z_cox_formula)) 
  res = coxph(formula_str, data = data_aug_cc) 
  
  return(res)
}

#' Time-Dependent Augmented Complete Case Analysis (Robust) 
#'
#' @param data A list containing:
#' #'   - `df`: a data frame with one row per subject and the following columns:
#'
#'     - `id`: unique subject identifier
#'     - `Z_1`, `Z_2`, ...: individual-level covariates
#'     - `Y`: observed time (either failure or censoring) assuming full linkage
#'     - `Delta`: failure indicator (1 = event, 0 = censored) in the linked dataset
#'     - `linkage`: linkage indicator (1 = linked to followup observational study, 0 = not linked)
#'     - `censored_in_trial`: binary indicator (1 = censored during trial period, 0 = otherwise)
#'     - `failed_in_trial`: binary indicator (1 = failed during trial period, 0 = otherwise)
#'     - `Y_1`: observed time in the trial-only data, assuming no linkage
#'     
#'   - `tdata`: a data frame in counting process format, with one or more rows per subject.
#'     It must include:
#'       - `id`: unique subject identifier
#'       - `tstart`, `tstop`: interval start and stop times
#'       - `death`: event indicator (1 = event, 0 = censored)
#'       - Covariates including any time-dependent ones such as `Z_3` 
#'       - t_0 (optional): change point in time
#'
#' @param cox_covariates A character vector of column names in `tdata` to be included as covariates in the Cox model.
#'
#' @return A `coxph` model object fit to  the linked (complete-case) data, using the augmented complete case analysis,
#'  with robust variance via clustering by subject ID.

tdc_augmented_complete_case_analysis_robust <- function(data, cox_covariates) {
  # augmented complete case use patients that are linked or failed in trial.
  df = data$df
  linked = which(df$linkage == 1) 
  failed_in_trial =  df$failed_in_trial  
  failed_in_trial_idx <- which(failed_in_trial == 1)
  linked_plus_failed_in_trial <- sort(unique(c(linked, failed_in_trial_idx))) 
  data_aug_cc = df[df$id %in% linked_plus_failed_in_trial, ]
  
  tdata = data$tdata
  tdata_aug_cc = tdata[tdata$id %in% linked_plus_failed_in_trial, ]
  
  # --- Cox model design
  Z_cox_formula <- paste(cox_covariates, collapse = " + ")
  formula_str <- as.formula(paste("Surv(tstart, tstop, death) ~", Z_cox_formula,
                                  ' + cluster(id)'))
  #res = coxph(Surv(tstart, tstop, death) ~ Z_1 + Z_2 + Z_3 + cluster(id), data=tdata_aug_cc)
  res = coxph(formula_str, data = tdata_aug_cc)
  return(res)
}

##############################################
##             NLAC analysis                ##
## For unlinked patients, censor them       ##
## at trial  (i.e., administrative cutoff)  ##
##############################################
 
#' Time-Independent Not Linked as Censored (NLAC) Analysis 
#'
#'@param data A list containing:

#'   - `df`: a data frame with one row per subject and the following columns:
#'
#'     - `id`: unique subject identifier
#'     - `Z_1`, `Z_2`, ...: individual-level covariates
#'     - `Y`: observed time (either failure or censoring) assuming full linkage
#'     - `Delta`: failure indicator (1 = event, 0 = censored) in the linked dataset
#'     - `linkage`: linkage indicator (1 = linked to followup observational study, 0 = not linked)
#'     - `censored_in_trial`: binary indicator (1 = censored during trial period, 0 = otherwise)
#'     - `failed_in_trial`: binary indicator (1 = failed during trial period, 0 = otherwise)
#'     - `Y_1`: observed time in the trial-only data, assuming no linkage
#'     
#' @param cox_covariates A character vector of column names in `df` to be included as covariates in the Cox model.
#' @return A `coxph` model object fit to the linked data, using the NLAC analysis.
#' 
ti_nlac_analysis <- function(data, cox_covariates) {
  df = data$df
  n = nrow(df)
  linked = which(df$linkage == 1)
  
  # Update Y and Delta for unlinked people
  df$Y[-linked] = df$Y_1[-linked]
  df$Delta[-linked] = df$failed_in_trial[-linked]
  
  # --- Cox model design
  Z_cox_formula <- paste(cox_covariates, collapse = " + ")
  formula_str <- as.formula(paste("Surv(Y, Delta) ~", Z_cox_formula)) 
  res = coxph(formula_str, data = df)
  
  #  res = coxph(Surv(Y, Delta) ~ Z_1 + Z_2 + Z_3, data = data)
  
  return(res)
}  

#' Time-Dependent Not Linked as Censored (NLAC) Analysis (Robust) 
#'
#' @param data A list containing:
#' #'   - `df`: a data frame with one row per subject and the following columns:
#'
#'     - `id`: unique subject identifier
#'     - `Z_1`, `Z_2`, ...: individual-level covariates
#'     - `Y`: observed time (either failure or censoring) assuming full linkage
#'     - `Delta`: failure indicator (1 = event, 0 = censored) in the linked dataset
#'     - `linkage`: linkage indicator (1 = linked to followup observational study, 0 = not linked)
#'     - `censored_in_trial`: binary indicator (1 = censored during trial period, 0 = otherwise)
#'     - `failed_in_trial`: binary indicator (1 = failed during trial period, 0 = otherwise)
#'     - `Y_1`: observed time in the trial-only data, assuming no linkage
#'     
#'   - `tdata`: a data frame in counting process format, with one or more rows per subject.
#'     It must include:
#'       - `id`: unique subject identifier
#'       - `tstart`, `tstop`: interval start and stop times
#'       - `death`: event indicator (1 = event, 0 = censored)
#'       - Covariates including any time-dependent ones such as `Z_3` 
#'       - t_0: change point in time
#'
#' @param cox_covariates A character vector of column names in `tdata` to be included as covariates in the Cox model.
#'
#' @return A `coxph` model object fit to  the linked (complete-case) data, using the NLAC analysis,
#'  with robust variance via clustering by subject ID.

tdc_nlac_analysis_robust <- function(data, cox_covariates) { 
  df = data$df
  tdata <- data$tdata
  n = nrow(df)
  linked = which(df$linkage == 1)
  
  # Update Y and Delta for unlinked people
  df$Y[-linked] = df$Y_1[-linked]
  df$Delta[-linked] = df$failed_in_trial[-linked]
  
  # create new tdata
  new_data = with(df, data.frame(id = 1:n,
                                   Z_1 = Z_1,
                                   Z_2 = Z_2, 
                                   Y = Y,
                                   Delta = Delta))
  tdata = tmerge(new_data, new_data, id=id, death=event(Y, Delta),
                 Z_3=tdc(df$t_0, df$Z_1))
  tdata$Z_3[is.na(tdata$Z_3)] <- 0
  
  # --- Cox model design
  Z_cox_formula <- paste(cox_covariates, collapse = " + ")
  formula_str <- as.formula(paste("Surv(tstart, tstop, death) ~", Z_cox_formula,
                                  ' + cluster(id)'))
  res = coxph(formula_str, data=tdata)
  return(res)
}

##############################################
##                                          ## 
##            IPLW analysis                 ## 
##                                          ## 
##############################################

#' Time-Independent Inverse Probability Weighted (IPW) Cox Analysis 
#' 
#'@param data A list containing:

#'   - `df`: a data frame with one row per subject and the following columns:
#'
#'     - `id`: unique subject identifier
#'     - `Z_1`, `Z_2`, ...: individual-level covariates
#'     - `Y`: observed time (either failure or censoring) assuming full linkage
#'     - `Delta`: failure indicator (1 = event, 0 = censored) in the linked dataset
#'     - `linkage`: linkage indicator (1 = linked to followup observational study, 0 = not linked)
#'     - `censored_in_trial`: binary indicator (1 = censored during trial period, 0 = otherwise)
#'     - `failed_in_trial`: binary indicator (1 = failed during trial period, 0 = otherwise)
#'     - `Y_1`: observed time in the trial-only data, assuming no linkage
#'     
#' @param cox_covariates A character vector of column names in `df` to be included as covariates in the Cox model.
#' @param link_covariates A character vector of column names in `df` to be included as covariates in the linkage model.
#' @return A list with:
#'   - `beta_hat`: estimated regression coefficients
#'   - `beta_hat_se`: robust standard errors
#'   - `beta_var`: sandwich variance estimator
#' 
ti_ipw_analysis <- function(data,
                            cox_covariates,
                            link_covariates) {
  expit <- function(x) 1 / (1 + exp(-x))
  
  # Subsets and design matrices
  df = data$df
  linkage <- df$linkage 
  Y <- df$Y
  Delta <- data$Delta 
  censored_in_trial <- df$censored_in_trial
  failed_in_trial <- df$failed_in_trial
  n <- length(Y)
  
  # --- Linkage model: P(linkage | Z) among censored in trial
  
  linked_and_censored_in_trial <- which(censored_in_trial == 1 & linkage == 1)
  
  formula_str <- paste("linkage ~", paste(link_covariates, collapse = " + ")) 
  gamma_glm <- glm(as.formula(formula_str),
                   data = df,
                   family = "binomial", subset = (censored_in_trial == 1))
  gamma_hat <- gamma_glm$coefficients
  Sigma_gamma_inv <- summary(gamma_glm)$cov.unscaled * n
  
  model_matrix_formula = paste("~ ", paste(link_covariates, collapse = " + "), sep="")
  Z_link_plus_intercept = model.matrix(as.formula(model_matrix_formula),
                                       data = df)
  estimated_prob <- expit(Z_link_plus_intercept[linked_and_censored_in_trial, ] %*% gamma_hat)[, 1]
  estimated_prob_full <- expit(Z_link_plus_intercept %*% gamma_hat)[, 1]
  
  # --- Set IP weights
  ip_weights <- rep(1, n)
  ip_weights[linked_and_censored_in_trial] <- 1 / estimated_prob
  
  # --- Cox model design 
  Z_cox_formula <- paste(cox_covariates, collapse = " + ")
  linked <- which(linkage == 1)
  failed_in_trial_idx <- which(failed_in_trial == 1)
  linked_plus_failed_in_trial <- sort(unique(c(linked, failed_in_trial_idx)))
  
  # --- Fit Cox model (non-TDC)
  data_temp <- cbind(data$df,ip_weights)
  data_linked_plus_failed_in_trial <- data_temp[linked_plus_failed_in_trial, ] 
  formula_str <- as.formula(paste("Surv(Y, Delta) ~", Z_cox_formula)) 
  res <- coxph(as.formula(formula_str),
               weights = ip_weights,
               data = data_linked_plus_failed_in_trial)
  
  beta_hat <- res$coefficients
  
  # --- Score residuals for influence function
  d = length(beta_hat)
  U_i_beta <- matrix(0, nrow = n, ncol = d)
  temp_score <- resid(res, type = "score",
                      collapse = data_temp$id[linked_plus_failed_in_trial],
                      weighted = FALSE)
  U_i_beta[linked_plus_failed_in_trial, ] <- temp_score
  
  Q_e_U_i_beta <- t(Z_link_plus_intercept[linked_and_censored_in_trial, ] *
                      (1 - estimated_prob) / estimated_prob) %*%
    U_i_beta[linked_and_censored_in_trial, ] / n
   
  # --- Sandwich variance 
  temp <- U_i_beta * ip_weights
  temp_1 <- Z_link_plus_intercept %*% Sigma_gamma_inv %*% Q_e_U_i_beta *
    censored_in_trial * (linkage - estimated_prob_full)
  temp_2 <- temp - temp_1
  Sigma <- t(temp_2) %*% temp_2 / n 
  
  Sigma_0_inv <- res$naive.var * n
  beta_var <- Sigma_0_inv %*% Sigma %*% Sigma_0_inv
  beta_hat_se <- sqrt(diag(beta_var) / n)
  
  return(list("beta_hat" = beta_hat, "beta_hat_se" = beta_hat_se, 'beta_var' = beta_var))
}


#' Time-Independent Inverse Probability Weighted (IPW) Cox Analysis (Robust) 
#'
#' @param data A list containing:
#' #'   - `df`: a data frame with one row per subject and the following columns:
#'
#'     - `id`: unique subject identifier
#'     - `Z_1`, `Z_2`, ...: individual-level covariates
#'     - `Y`: observed time (either failure or censoring) assuming full linkage
#'     - `Delta`: failure indicator (1 = event, 0 = censored) in the linked dataset
#'     - `linkage`: linkage indicator (1 = linked to followup observational study, 0 = not linked)
#'     - `censored_in_trial`: binary indicator (1 = censored during trial period, 0 = otherwise)
#'     - `failed_in_trial`: binary indicator (1 = failed during trial period, 0 = otherwise)
#'     - `Y_1`: observed time in the trial-only data, assuming no linkage
#'     
#'   - `tdata`: a data frame in counting process format, with one or more rows per subject.
#'     It must include:
#'       - `id`: unique subject identifier
#'       - `tstart`, `tstop`: interval start and stop times
#'       - `death`: event indicator (1 = event, 0 = censored)
#'       - Covariates including any time-dependent ones such as `Z_3` 
#'       - t_0 (optional): change point in time
#'
#' @param cox_covariates A character vector of column names in `tdata` to be included as covariates in the Cox model.
#' @param link_covariates A character vector of column names in `df` to be included as covariates in the linkage model.

#' @return A list with:
#'   - `beta_hat`: estimated regression coefficients
#'   - `beta_hat_se`: robust standard errors
#'   - `beta_var`: sandwich variance estimator

tdc_ipw_analysis <- function(data, cox_covariates, 
                             link_covariates) {
  # Subsets and design matrices
  df = data$df
  linkage <- df$linkage 
  Y <- df$Y
  Delta <- data$Delta 
  censored_in_trial <- df$censored_in_trial
  failed_in_trial <- df$failed_in_trial
  n <- length(Y)
  
  # --- Linkage model: P(linkage | Z) among censored in trial
  
  linked_and_censored_in_trial <- which(censored_in_trial == 1 & linkage == 1)
  
  formula_str <- paste("linkage ~", paste(link_covariates, collapse = " + ")) 
  gamma_glm <- glm(as.formula(formula_str),
                   data = df,
                   family = "binomial", subset = (censored_in_trial == 1))
  gamma_hat <- gamma_glm$coefficients
  Sigma_gamma_inv <- summary(gamma_glm)$cov.unscaled * n
  
  model_matrix_formula = paste("~ ", paste(link_covariates, collapse = " + "), sep="")
  Z_link_plus_intercept = model.matrix(as.formula(model_matrix_formula),
                                  data = df)
  estimated_prob <- expit(Z_link_plus_intercept[linked_and_censored_in_trial, ] %*% gamma_hat)[, 1]
  estimated_prob_full <- expit(Z_link_plus_intercept %*% gamma_hat)[, 1]
  
  # --- Set IP weights
  ip_weights <- rep(1, n)
  ip_weights[linked_and_censored_in_trial] <- 1 / estimated_prob
  
  # --- Cox model design
  Z_cox_formula <- paste(cox_covariates, collapse = " + ")
  
  linked <- which(linkage == 1)
  failed_in_trial_idx <- which(failed_in_trial == 1)
  linked_plus_failed_in_trial <- sort(unique(c(linked, failed_in_trial_idx)))
  
  tdata = data$tdata
  tdata_acc = tdata[tdata$id %in% linked_plus_failed_in_trial, ] 
  formula_str <- as.formula(paste("Surv(tstart, tstop, death) ~", Z_cox_formula))
  res <- coxph(formula_str, data = tdata_acc, weights = ip_weights[tdata_acc$id])
  beta_hat = res$coefficients
  features = names(beta_hat)
  
  # Score residuals 
  d = length(beta_hat)
  U_i_beta = matrix(0, nrow=n, ncol=d)
  temp = resid(res, type='score', collapse=tdata_acc$id, weighted=F)
  U_i_beta[linked_plus_failed_in_trial, ] = temp
  
  # --- Estimate cross-terms for influence function
  Q_e_U_i_beta <- t(Z_link_plus_intercept[linked_and_censored_in_trial, ] *
                      (1 - estimated_prob) / estimated_prob) %*%
    U_i_beta[linked_and_censored_in_trial, ] / n
  
  # --- Sandwich variance estimation 
  temp <- U_i_beta * ip_weights
  temp_1 <- Z_link_plus_intercept %*% Sigma_gamma_inv %*% Q_e_U_i_beta * 
    censored_in_trial * (linkage - estimated_prob_full)
  temp_2 <- temp - temp_1
  Sigma <- t(temp_2) %*% temp_2 / n  
 
  if (is.null(res$naive.var)) {
    Sigma_0_inv = res$var * n
  } else {
    Sigma_0_inv = res$naive.var * n
  }
   
  beta_var <- Sigma_0_inv %*% Sigma %*% Sigma_0_inv
  beta_hat <- res$coefficients
  beta_hat_se <- sqrt(diag(beta_var) / n)
  
  return(list("beta_hat"=beta_hat, "beta_hat_se"=beta_hat_se, 'beta_var' = beta_var))
}


##############################################
##                                          ## 
##           AIPLW analysis                 ## 
##                                          ## 
##############################################
# Functions that are necessary for aiplw estimators. 
get_risk_set <- function(
    time_of_interest,
    observed_time_vector
) {
  return(
    which(observed_time_vector >= time_of_interest)
  )
}

S_n_part <- function(t, beta, beta_observed, X, Y, weight, censored_in_trial_idx,
                     linked_plus_failed_in_trial, n, baseline_cum) {
  # Compute S_n_DR(beta, t) as 
  # S_n_DR(beta, t) = S_n_DR_1(beta, t) / S_n_DR_0(beta, t)
  # More specifically, we have 
  # S_n_DR_1(beta, t) = 1 / n \sum_{j=1}^n [hat w_j I(tilde T_j >= t) exp(X_j^T beta) + 
  # I(Q_j = 0) - hat w_j hat m_3(X_j, t)] X_j
  # S_n_DR_0(beta, t) = 1 / n \sum_{j=1}^n [hat w_j I(tilde T_j >= t) exp(X_j^T beta) + 
  # I(Q_j = 0) - hat w_j hat m_3(X_j, t)]
  # 
  # Args
  # t: the time to compute the S_n_DR(beta, t)
  # beta: beta parameter
  # beta_observed: beta for the cox model for observed time, used to estimate m_3. 
  # X: X matrix
  # Y: Y, only for linked_plus_failed_in_trial
  # weight: weight vector, full length
  # censored_in_trial_idx: index for patients that are censored in trial.
  # linked_plus_failed_in_trial: index for patients that are linked or failed in trial.
  # n: number of samples
  # baseline_cum: baseline cumulative hazards for the observed time at time t. need this to 
  # compute the m_3(X_j, t).
  # m_3(X_j, t) = P(T_tilde >= t | X_j, Q = 0, L = 1) exp(X_j^T beta)
  #             = exp(-Lambda_0(t) exp(X_j^T beta_observed)) exp(X_j^T beta)
  risk_set = get_risk_set(t, Y)
  X_risk_set = X[linked_plus_failed_in_trial, ][risk_set, ]
  temp_result = exp(X_risk_set %*% beta) * weight[linked_plus_failed_in_trial][risk_set]
  t1 = t(X_risk_set) %*% temp_result
  t1_denom = sum(temp_result)
  
  # The following part seems to be affecting the estimating of U_n_AIPW
  t2_weight = 1 - weight[censored_in_trial_idx]
  X_censored_in_trial = X[censored_in_trial_idx, ]
  exp_X_beta = exp(X_censored_in_trial %*% beta)
  exp_X_beta_observed = exp(X_censored_in_trial %*% beta_observed)
  m_3 = exp(-baseline_cum * exp_X_beta_observed) * exp_X_beta
  temp_result = m_3 * t2_weight
  t2 = t(X_censored_in_trial) %*% temp_result
  t2_denom = sum(temp_result)
  return((t1 + t2) / (t1_denom + t2_denom))
}

S_n_covariance <- function(
    t, beta, beta_observed, X, Y, weight, censored_in_trial_idx,
    linked_plus_failed_in_trial, n, baseline_cum
) {
  # Compute [grad_beta S_n_DR(beta, t)^1] * S_n_DR(beta, t)^0 - S_n_DR(beta, t)^1 * [grad_beta S_n_DR(beta, t)^0]
  # divided by [S_n_DR(beta, t)^0]^2
  # This is used to compute the gradient of U_n_AIPW.
  
  # Args
  # t: the time to compute the S_n_DR(beta, t)
  # beta: beta parameter
  # beta_observed: beta for the cox model for observed time, used to estimate m_3. 
  # X: X matrix
  # Y: Y, only for linked_plus_failed_in_trial
  # weight: weight vector, full length
  # censored_in_trial_idx: index for patients that are censored in trial.
  # linked_plus_failed_in_trial: index for patients that are linked or failed in trial.
  # n: number of samples
  # baseline_cum: baseline cumulative hazards for the observed time at time t. need this to 
  # compute the m_3(X_j, t).
  # m_3(X_j, t) = P(T_tilde >= t | X_j, Q = 0, L = 1) exp(X_j^T beta)
  #             = exp(-Lambda_0(t) exp(X_j^T beta)) exp(X_j^T beta)
  # index_set: subindex for linked_plus_failed_in_trial that are linked_and_censored_in_trial. 
  
  # First, we compute the [grad_beta S_n_DR(beta, t)^1]. 
  
  risk_set = get_risk_set(t, Y)
  X_risk_set = X[linked_plus_failed_in_trial, ][risk_set, ]
  temp_result = exp(X_risk_set %*% beta) * weight[linked_plus_failed_in_trial][risk_set]
  S_0_part_1 = sum(temp_result) # used to compute S_0 later.
  S_1_part_1 = t(X_risk_set) %*% temp_result # used to compute S_1 later.
  
  X_risk_set_scale = TRA(t(X_risk_set), sqrt(temp_result), '*')
  S_1_grad_part_1 = X_risk_set_scale %*% t(X_risk_set_scale)
  
  t2_weight = 1 - weight[censored_in_trial_idx]
  X_censored_in_trial = X[censored_in_trial_idx, ]
  exp_X_beta = exp(X_censored_in_trial %*% beta) 
  exp_X_beta_observed = exp(X_censored_in_trial %*% beta_observed)
  m_3 = exp(-baseline_cum * exp_X_beta_observed) * exp_X_beta
  m_3_t2_weight = m_3 * t2_weight
  S_0_part_2 = sum(m_3_t2_weight) # used to compute S_0 later
  S_1_part_2 = t(X_censored_in_trial) %*% m_3_t2_weight # used to compute S_1 later.
  
  X_temp_result = TRA(t(X_censored_in_trial), m_3, '*')
  X_temp_result_2 = TRA(t(X_censored_in_trial), t2_weight, '*')
  S_1_grad_part_2 = X_temp_result %*% t(X_temp_result_2)
  
  S_1_grad = (S_1_grad_part_1 + S_1_grad_part_2) / n
  S_1_grad = (S_1_grad + t(S_1_grad)) / 2
  
  # Next, we compute the S_n_DR(beta, t)^0.
  S_0 = (S_0_part_1 + S_0_part_2) / n
  
  # Further, we compute S_n_DR(beta, t)^1.
  S_1 = (S_1_part_1 + S_1_part_2) / n
  
  # Now compute [grad_beta S_n_DR(beta, t)^0], which is actually the same thing as S_n_DR(beta, t)^1. 
  # Thus, we have that for the first part of grad_beta U_n_AIPW
  result = S_1_grad / S_0 - S_1 %*% t(S_1) / (S_0^2)
  
  return(result)
}

covariance_compute_aipw <- function(
    beta, beta_observed, X, Y, delta, weight, n, censored_in_trial_idx,
    linked_plus_failed_in_trial, linked_and_censored_in_trial, 
    baseline_ch, estimated_m1, Y_response, Delta_3,
    Z_plus_intercept, index_set
) {
  # Args
  # beta: beta parameter
  # beta_observed: beta for the cox model for observed time, used to estimate m_3. 
  # X: X matrix, should contain all the predictors.
  # Y: Y, only for linked_plus_failed_in_trial.
  # delta: delta, also only for linked_plus_failed_in_trial. 
  # weight: weight vector, full length.
  # n: full number of observations.
  # censored_in_trial_idx: stores the indices for censored_in_trial.
  # linked_plus_failed_in_trial: stores the indices for linked_plus_failed_in_trial.
  # baseline_ch: baseline cumulative hazards for observed time T tilde. 
  # estimated_m1: estimated values for E[Delta | X, Q = 0, L = 1], only for censored_in_trial.
  # Y_response: Y, only for linked_and_censored_in_trial.
  # Delta_3: delta, only for linked_and_censored_in_trial.
  # Z_plus_intercept: intercept + X matrix. 
  # index_set: subindex for linked_plus_failed_in_trial that are linked_and_censored_in_trial. 
  # Used to speed up the computation for gradient. 
  
  t_data = cbind(Y, baseline_ch[linked_plus_failed_in_trial])
  
  weight_delta = matrix(weight[linked_plus_failed_in_trial] * delta, ncol=1)
  covariance_first_part <- apply(
    t_data, 1,
    function(df) {
      ts <- df[1]
      baseline_cum <- df[2]
      s_n_cov <- S_n_covariance(
        ts, beta, beta_observed, X, Y, weight, censored_in_trial_idx, linked_plus_failed_in_trial,
        n, baseline_cum
      )
      return(s_n_cov)
    }
  )
  first_part = covariance_first_part %*% weight_delta / n
  first_part = matrix(first_part, nrow=3)
  
  # Now look at the second part 
  second_part_weight = 1 - weight[censored_in_trial_idx]
  m2_grad_response = covariance_first_part[, index_set]
  
  # Use linear regression to estimate the three components
  m2_grad_reg_1 = lm(m2_grad_response[1, ] * Delta_3 ~ 
                       X[linked_and_censored_in_trial, 1] + 
                       X[linked_and_censored_in_trial, 2] + 
                       X[linked_and_censored_in_trial, 3])
  
  m2_grad_reg_2 = lm(m2_grad_response[2, ] * Delta_3 ~ 
                       X[linked_and_censored_in_trial, 1] + 
                       X[linked_and_censored_in_trial, 2] + 
                       X[linked_and_censored_in_trial, 3])
  
  m2_grad_reg_3 = lm(m2_grad_response[3, ] * Delta_3 ~ 
                       X[linked_and_censored_in_trial, 1] + 
                       X[linked_and_censored_in_trial, 2] + 
                       X[linked_and_censored_in_trial, 3])
  
  m2_grad_reg_5 = lm(m2_grad_response[5, ] * Delta_3 ~ 
                       X[linked_and_censored_in_trial, 1] + 
                       X[linked_and_censored_in_trial, 2] + 
                       X[linked_and_censored_in_trial, 3])
  
  m2_grad_reg_6 = lm(m2_grad_response[6, ] * Delta_3 ~ 
                       X[linked_and_censored_in_trial, 1] + 
                       X[linked_and_censored_in_trial, 2] + 
                       X[linked_and_censored_in_trial, 3])
  
  m2_grad_reg_9 = lm(m2_grad_response[9, ] * Delta_3 ~ 
                       X[linked_and_censored_in_trial, 1] + 
                       X[linked_and_censored_in_trial, 2] + 
                       X[linked_and_censored_in_trial, 3])
  
  m2_comp_1 = Z_plus_intercept[censored_in_trial_idx, ] %*% m2_grad_reg_1$coefficients
  m2_comp_2 = Z_plus_intercept[censored_in_trial_idx, ] %*% m2_grad_reg_2$coefficients
  m2_comp_3 = Z_plus_intercept[censored_in_trial_idx, ] %*% m2_grad_reg_3$coefficients
  m2_comp_5 = Z_plus_intercept[censored_in_trial_idx, ] %*% m2_grad_reg_5$coefficients
  m2_comp_6 = Z_plus_intercept[censored_in_trial_idx, ] %*% m2_grad_reg_6$coefficients
  m2_comp_9 = Z_plus_intercept[censored_in_trial_idx, ] %*% m2_grad_reg_9$coefficients
  
  second_part_grad = t(
    cbind(
      m2_comp_1, m2_comp_2, m2_comp_3,
      m2_comp_2, m2_comp_5, m2_comp_6,
      m2_comp_3, m2_comp_6, m2_comp_9
    )
  ) %*% second_part_weight
  second_part = -matrix(second_part_grad, nrow=3) / n
  cov_est = first_part - second_part
  cov_est = solve(cov_est)
  return(cov_est)  
}

gradient_compute_aipw <- function(beta, beta_observed, X, Y, delta, weight, n, censored_in_trial_idx,
                                  linked_plus_failed_in_trial, linked_and_censored_in_trial, 
                                  baseline_ch, estimated_m1, Y_response, Delta_3,
                                  Z_plus_intercept, index_set) {
  # Args
  # beta: beta parameter
  # beta_observed: beta parameter of cox model for observed time. 
  # X: X matrix, should contain all the predictors.
  # Y: Y, only for linked_plus_failed_in_trial.
  # delta: delta, also only for linked_plus_failed_in_trial. 
  # weight: weight vector, full length.
  # n: full number of observations.
  # censored_in_trial_idx: stores the indices for censored_in_trial.
  # linked_plus_failed_in_trial: stores the indices for linked_plus_failed_in_trial.
  # baseline_ch: baseline cumulative hazards for observed time T tilde. 
  # estimated_m1: estimated values for E[Delta | X, Q = 0, L = 1], only for censored_in_trial.
  # Y_response: Y, only for linked_and_censored_in_trial.
  # Delta_3: delta, only for linked_and_censored_in_trial.
  # Z_plus_intercept: intercept + X matrix. 
  # index_set: subindex for linked_plus_failed_in_trial that are linked_and_censored_in_trial. 
  # Used to speed up the computation for gradient. 
  
  t_data = cbind(Y, baseline_ch[linked_plus_failed_in_trial])
  
  weight_delta = matrix(weight[linked_plus_failed_in_trial] * delta, ncol=1)
  
  gradient_first_part <- apply(t_data, 1,
                               function(df) {
                                 ts <- df[1]
                                 baseline_cum <- df[2]
                                 s_n <- S_n_part(ts, beta, 
                                                 beta_observed,
                                                 X, Y, weight, censored_in_trial_idx, 
                                                 linked_plus_failed_in_trial,
                                                 n, baseline_cum)
                                 return(s_n)
                               })
  # first_part = gradient_first_part %*% weight_delta / n
  first_part = (t(X[linked_plus_failed_in_trial, ]) - gradient_first_part) %*% weight_delta / n
  
  second_part_weight = 1 - weight[censored_in_trial_idx]
  second_part_weight_1 = second_part_weight * estimated_m1
  m_0_first_part = t(X[censored_in_trial_idx, ]) %*% second_part_weight_1
   
  m2_response = gradient_first_part[, index_set]
  
  # Use linear regression to estimate the three components
  m2_reg_1 = lm(m2_response[1, ] * Delta_3 ~
                  X[linked_and_censored_in_trial, 1] +
                  X[linked_and_censored_in_trial, 2] +
                  X[linked_and_censored_in_trial, 3])
  
  m2_reg_2 = lm(m2_response[2, ] * Delta_3 ~
                  X[linked_and_censored_in_trial, 1] +
                  X[linked_and_censored_in_trial, 2] +
                  X[linked_and_censored_in_trial, 3])
  
  m2_reg_3 = lm(m2_response[3, ] * Delta_3 ~
                  X[linked_and_censored_in_trial, 1] +
                  X[linked_and_censored_in_trial, 2] +
                  X[linked_and_censored_in_trial, 3])
  
  m2_comp_1 = Z_plus_intercept[censored_in_trial_idx, ] %*% m2_reg_1$coefficients
  m2_comp_2 = Z_plus_intercept[censored_in_trial_idx, ] %*% m2_reg_2$coefficients
  m2_comp_3 = Z_plus_intercept[censored_in_trial_idx, ] %*% m2_reg_3$coefficients
  m_0_second_part = t(cbind(m2_comp_1, m2_comp_2, m2_comp_3)) %*% second_part_weight
  second_part = (m_0_first_part - m_0_second_part) / n
  
  return(first_part + second_part)
}

aipw_estimate <- function(beta, data) {
  # AIPLW estimator for beta
  # Here we try to solve estimating equation 17 from paper 
  # Args: data
  df = data$df
  linkage = df$linkage
  n = length(linkage)
  linked = which(df$linkage == 1)
  
  Y = df$Y
  Delta = df$Delta
  
  Z_1 = df$Z_1
  Z_2 = df$Z_2
  Z_3 = df$Z_3
   
  censored_in_trial = df$censored_in_trial
  failed_in_trial = df$failed_in_trial
  
  censored_in_trial_idx = which(censored_in_trial == 1)
  # Estimate the linkage probability and weights.
  # We use logistic regression model to estimate the linkage probability. 
  gamma_glm = glm(linkage[censored_in_trial_idx] ~ Z_1[censored_in_trial_idx] + 
                    Z_2[censored_in_trial_idx] + Z_3[censored_in_trial_idx],
                  family="binomial")
  gamma_hat = gamma_glm$coefficients
  
  Z = cbind(Z_1, Z_2, Z_3)
  Z_plus_intercept = cbind(rep(1, n), Z)
  linked_and_censored_in_trial = which(censored_in_trial == 1
                                       & linkage == 1)
  estimated_prob = expit(Z_plus_intercept[linked_and_censored_in_trial, ] %*% gamma_hat)[, 1]
  
  ip_weights = rep(0, n)
  failed_in_trial_idx = which(failed_in_trial == 1)
  ip_weights[failed_in_trial_idx] = 1
  ip_weights[linked_and_censored_in_trial] = 1 / estimated_prob
  
  # Estimate m_3(X, t) = P(Y >= t | X, Q = 0, L = 1) exp(X^T beta).
  # We use cox regression model to estimate the survival probability for Y.
  # and Y is the observed time. 
  # With cox model, we have P(Y >= t | X, Q = 0, L = 1) = 
  # exp(-Lambda_0(t) exp(X^T beta_observed)).
  data_combine = cbind(
    Y[linked_and_censored_in_trial],
    rep(1, length(linked_and_censored_in_trial)), # we are computing survival probability for observed time, thus we can treat it as no-censoring for everyone. 
    Z_1[linked_and_censored_in_trial],
    Z_2[linked_and_censored_in_trial],
    Z_3[linked_and_censored_in_trial],
    Delta[linked_and_censored_in_trial]
  )
  
  data_combine = as.data.frame(data_combine)
  colnames(data_combine) = c("Y", "Delta2", "Z1", "Z2", "Z3", "Delta")
  surv_observed = coxph(Surv(Y, Delta2) ~ Z1 + Z2 + Z3, data=data_combine)
  beta_observed = surv_observed$coefficients
  
  H0 <- basehaz(surv_observed, centered=FALSE)  # centered False means that we are setting all covariates values to zero. 
  baseline_ch_fun = stepfun(H0$time, c(0, H0$hazard))
  baseline_ch = baseline_ch_fun(Y)
  
  # now estimate m_1(x) = E[Delta | X, Q = 0, L = 1] X, we use logistic regression to estimate the 
  # E[Delta | X, Q = 0, L = 1] part. 
  m1_glm = glm(Delta ~ Z1 + Z2 + Z3, data=data_combine, family="binomial")
  m1_glm_hat = m1_glm$coefficients
  estimated_m1 = expit(Z_plus_intercept[censored_in_trial_idx, ] %*% m1_glm_hat)[, 1]
  
  linked_plus_failed_in_trial = sort(unique(c(linked, failed_in_trial_idx)))
  Y_response = Y[linked_and_censored_in_trial]
  Delta_3 = Delta[linked_and_censored_in_trial]
  
  # Try to get the index of elements that are linked and also censored in trial
  # within linked_plus_failed_in_trial index set. 
  #Q = as.numeric(data$C_PCPT > data$T_failure)
  linkage_restricted = linkage[linked_plus_failed_in_trial]
  Q_restricted = failed_in_trial[linked_plus_failed_in_trial]
  index_set = which(linkage_restricted == 1 & Q_restricted == 0)
  
  n_iters = 5
  store_beta = matrix(0, nrow=n_iters, ncol=3)
  
  start_time = Sys.time()
  for (i in 1:n_iters) {
    # if (i == 1) {
    #   cat("beta has value", beta, "\n")
    # }
    store_beta[i, ] = beta
    gradient <- gradient_compute_aipw(
      beta, beta_observed, Z, Y[linked_plus_failed_in_trial], 
      Delta[linked_plus_failed_in_trial], ip_weights, 
      n, censored_in_trial_idx,
      linked_plus_failed_in_trial, linked_and_censored_in_trial,
      baseline_ch, estimated_m1,
      Y_response, Delta_3, Z_plus_intercept, index_set
    )
    cov_est <- covariance_compute_aipw(
      beta, beta_observed, Z, Y[linked_plus_failed_in_trial], 
      Delta[linked_plus_failed_in_trial], ip_weights, 
      n, censored_in_trial_idx,
      linked_plus_failed_in_trial, linked_and_censored_in_trial,
      baseline_ch, estimated_m1,
      Y_response, Delta_3, Z_plus_intercept, index_set
    )
    gradient_norm = norm(gradient, type="2")
    # cat("the gradient has value", norm(gradient, type="2"), "\n\n")
    
    if (gradient_norm < 1e-4) {
      break
    }
    beta = beta + cov_est %*% gradient
    # cat("beta has value", beta, "\n")
  }
  return(beta)
}

aipw_bootstrap <- function(data, n_BT) {
  # Here we use bootstrap to estimate the variance of our aipw estimates.
  df = data$df
  n = length(df$Y)
  beta_initial = rep(0, 3)
  beta_estimate = aipw_estimate(beta_initial, data)
  beta_bt_estimate = matrix(0, nrow=n_BT, ncol=3)
  
  for (i in 1:n_BT) {
    if (i %% 10 == 0) {
      cat("now is the ", i, " iteration\n")
    }
    bootstrap_index = sample(n, size=n, replace=T)
    data_BT = list()
    data_BT$df = data$df[bootstrap_index,] 
    beta_bt_estimate[i, ] = aipw_estimate(beta_estimate, data_BT)
  }
  return(list(
    "beta_estimate"=beta_estimate,
    "beta_bt_estimate"=beta_bt_estimate
  )) 
}

all_simulation <- function(n, seed, method_num, 
                           setting = 1, tdcAnalysis = FALSE, 
                           robust = TRUE, n_BT = 100) {
  set.seed(seed)
  if (setting == 1) {
    data = tdc_data_generation_setting_1(n, seed)
    beta_1 = -log(4)
    beta_2 = log(1.5)
    beta_3 = 0.5
    robust = TRUE 
    tdc = TRUE
  } else if (setting == 2){
    data = tdc_data_generation_setting_2(n, seed)
    beta_1 = -log(4)
    beta_2 = log(1.5)
    beta_3 = 0.5
    robust = TRUE 
    tdc = TRUE
  } else if (setting == 3){
    data = tdc_data_generation_setting_3(n, seed)
    beta_1 = -log(4)
    beta_2 = log(1.5)
    beta_3 = 0.5
    robust = TRUE 
    tdc = TRUE
  } else if (setting == 4){
    data = tdc_data_generation_setting_4(n, seed)
    beta_1 = -log(4)
    beta_2 = log(1.5)
    beta_3 = 0.5
    robust = TRUE 
    tdc = TRUE
  } else if (setting == 5){
    data = tdc_data_generation_setting_5(n, seed)
    beta_1 = -log(4)
    beta_2 = log(1.5)
    beta_3 = 0.5
    robust = TRUE 
    tdc = TRUE
  }   
  if (setting == 6){
    data = ti_data_generation_setting_6(n, seed)
    beta_1=-log(2)
    beta_2=log(2)  
    beta_3 = 0.2
  } else if (setting == 7){
    data = ti_data_generation_setting_7(n, seed)
    beta_1=-log(2)
    beta_2=log(2)  
    beta_3 = 0.2
  }
  
  if (setting %in% c(1, 2, 3, 4, 5) & tdcAnalysis == TRUE){
    if (method_num == 1) {
      res_fh = tdc_full_cohort_analysis_robust(data, 
                                               cox_covariates = 
                                                 c("Z_1", "Z_2", "Z_3"))
      beta_hat = res_fh$coefficients
      beta_hat_se = sqrt(diag(res_fh$var))
    } else if (method_num == 2) {
      res_cc = tdc_complete_case_analysis_robust(data, 
                                                 cox_covariates = 
                                                   c("Z_1", "Z_2", "Z_3"))
      beta_hat = res_cc$coefficients
      beta_hat_se = sqrt(diag(res_cc$var))
    } else if (method_num == 3) {
      res_acc = tdc_augmented_complete_case_analysis_robust(data, 
                                                            cox_covariates = 
                                                              c("Z_1", "Z_2", "Z_3"))
      beta_hat = res_acc$coefficients
      beta_hat_se = sqrt(diag(res_acc$var))
    } else if (method_num == 4) {
      # For people that are not linked and Y > tau, set their Y to be tau and Delta to be 0.
      res_nlac = tdc_nlac_analysis_robust(data, cox_covariates = 
                                            c("Z_1", "Z_2", "Z_3")) 
      beta_hat = res_nlac$coefficients
      beta_hat_se = sqrt(diag(res_nlac$var))
    } else if (method_num == 5) { 
      res_ipw = tdc_ipw_analysis(data, 
                                 cox_covariates = c("Z_1", "Z_2", "Z_3"),
                                 link_covariates = c("Z_1", "Z_2"))
      beta_hat = res_ipw$beta_hat
      beta_hat_se = res_ipw$beta_hat_se
    }
  }
  
  if (setting %in% c(1, 2, 3, 5) & tdcAnalysis == FALSE){
    if (method_num == 1) { 
      beta_1 = -1.1
      res_fh = ti_full_cohort_analysis(data, cox_covariates = 
                                         c("Z_1", "Z_2"))
      beta_hat = res_fh$coefficients
      beta_hat_se = sqrt(diag(res_fh$var))
    } else if (method_num == 2) {
      beta_1 = -1.1
      res_cc = ti_complete_case_analysis(data, cox_covariates = 
                                           c("Z_1", "Z_2"))
      beta_hat = res_cc$coefficients
      beta_hat_se = sqrt(diag(res_cc$var))
    } else if (method_num == 3) {
      beta_1 = -1.1
      res_acc = ti_augmented_complete_case_analysis(data, cox_covariates = 
                                                      c("Z_1", "Z_2"))
      beta_hat = res_acc$coefficients
      beta_hat_se = sqrt(diag(res_acc$var))
    } else if (method_num == 4) {
      # For people that are not linked and Y > tau, set their Y to be tau and Delta to be 0.
      # Works.
      beta_1 = -1.1
      res_nlac = ti_nlac_analysis(data, cox_covariates = 
                                    c("Z_1", "Z_2"))
      beta_hat = res_nlac$coefficients
      beta_hat_se = sqrt(diag(res_nlac$var))
    } else if (method_num == 5) {
      beta_1 = -1.1
      res_ipw = ti_ipw_analysis(data, 
                                cox_covariates = c("Z_1", "Z_2"),
                                link_covariates = c("Z_1", "Z_2"))
      beta_hat = res_ipw$beta_hat
      beta_hat_se = res_ipw$beta_hat_se
    }
  }
  
  if (setting %in% c(6, 7)){
    if (method_num == 1) {
      res_fh = ti_full_cohort_analysis(data, cox_covariates = 
                                         c("Z_1", "Z_2", "Z_3"))
      beta_hat = res_fh$coefficients
      beta_hat_se = sqrt(diag(res_fh$var))
    } else if (method_num == 2) {
      res_cc = ti_complete_case_analysis(data, cox_covariates = 
                                           c("Z_1", "Z_2", "Z_3"))
      beta_hat = res_cc$coefficients
      beta_hat_se = sqrt(diag(res_cc$var))
    } else if (method_num == 3) {
      res_acc = ti_augmented_complete_case_analysis(data, cox_covariates = 
                                                      c("Z_1", "Z_2", "Z_3"))
      beta_hat = res_acc$coefficients
      beta_hat_se = sqrt(diag(res_acc$var))
    } else if (method_num == 4) {
      # For people that are not linked and Y > tau, set their Y to be tau and Delta to be 0.
      # Works.
      res_nlac = ti_nlac_analysis(data, cox_covariates = 
                                    c("Z_1", "Z_2", "Z_3"))
      beta_hat = res_nlac$coefficients
      beta_hat_se = sqrt(diag(res_nlac$var))
    } else if (method_num == 5) {
      res_ipw = ti_ipw_analysis(data, 
                                cox_covariates = c("Z_1", "Z_2", "Z_3"),
                                link_covariates = c("Z_1", "Z_2") )
      beta_hat = res_ipw$beta_hat
      beta_hat_se = res_ipw$beta_hat_se
    }  else if (method_num == 6) {
      start_time = Sys.time()
      result_aipw = aipw_bootstrap(data, n_BT)
      end_time = Sys.time()
      cat("time spend", end_time - start_time, "\n")
      beta_hat = result_aipw$beta_estimate
      beta_bt_estimate = result_aipw$beta_bt_estimate
      beta_hat_se = apply(beta_bt_estimate, 2, sd)
    }
  }
  
  covered_beta_1 = beta_hat[1] - beta_hat_se[1] * 1.96 <= beta_1 & beta_1 <= beta_hat[1] + beta_hat_se[1] * 1.96
  covered_beta_2 = beta_hat[2] - beta_hat_se[2] * 1.96 <= beta_2 & beta_2 <= beta_hat[2] + beta_hat_se[2] * 1.96
  covered_beta_3 = beta_hat[3] - beta_hat_se[3] * 1.96 <= beta_3 & beta_3 <= beta_hat[3] + beta_hat_se[3] * 1.96
  
  return(list("beta" = beta_hat,  "covered_beta_1" = covered_beta_1, "covered_beta_2" = covered_beta_2,
              "covered_beta_3" = covered_beta_3, "beta_hat_se" = beta_hat_se))
}

