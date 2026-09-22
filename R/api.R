#' Simulate linked trial and observational follow-up data
#'
#' Generate one dataset from a simulation design used in Cheng et al. (2026).
#'
#' @param n Positive integer sample size.
#' @param setting Integer from 1 through 7 selecting a paper simulation design.
#'   Settings 1--5 use a time-dependent treatment effect; settings 6--7 use
#'   time-independent covariates.
#' @param seed Integer random seed.
#'
#' @return A list with `df`, the subject-level data, and (for settings 1--5)
#'   `tdata`, a counting-process representation.
#' @export
simulate_linkage_data <- function(n, setting = 1L, seed = 1L) {
  .check_positive_integer(n, "n")
  .check_positive_integer(seed, "seed")
  setting <- as.integer(setting)
  if (length(setting) != 1L || is.na(setting) || !setting %in% 1:7) {
    stop("`setting` must be one integer from 1 through 7.", call. = FALSE)
  }

  generators <- list(
    tdc_data_generation_setting_1,
    tdc_data_generation_setting_2,
    tdc_data_generation_setting_3,
    tdc_data_generation_setting_4,
    tdc_data_generation_setting_5,
    ti_data_generation_setting_6,
    ti_data_generation_setting_7
  )
  generators[[setting]](as.integer(n), as.integer(seed))
}

#' Fit a Cox model under a linkage-analysis strategy
#'
#' Fit the oracle, complete-case (CC), augmented complete-case (CC+),
#' non-linked-as-censored (NLAC), or inverse-probability-of-linkage weighted
#' (IPLW) analysis described by Cheng et al. (2026).
#'
#' @param data A list containing a subject-level `df` data frame and, for a
#'   time-dependent analysis, a counting-process `tdata` data frame.
#' @param method One of `"oracle"`, `"complete_case"`,
#'   `"augmented_complete_case"`, `"nlac"`, or `"iplw"`.
#' @param cox_covariates Character vector naming Cox-model covariates.
#' @param linkage_covariates Character vector naming linkage-model covariates;
#'   required for `method = "iplw"`.
#' @param time_dependent Whether to fit a counting-process Cox model.
#'
#' @return An object of class `longterm_fit` with coefficients, standard
#'   errors, covariance matrix, and the underlying model result.
#' @export
fit_linkage_cox <- function(
    data,
    method = c(
      "iplw", "nlac", "complete_case", "augmented_complete_case", "oracle"
    ),
    cox_covariates,
    linkage_covariates = NULL,
    time_dependent = FALSE) {
  method <- match.arg(method)
  .validate_analysis_data(data, time_dependent)
  .check_covariates(cox_covariates, "cox_covariates")

  if (identical(method, "iplw")) {
    .check_covariates(linkage_covariates, "linkage_covariates")
  }

  result <- if (isTRUE(time_dependent)) {
    switch(
      method,
      oracle = tdc_full_cohort_analysis_robust(data, cox_covariates),
      complete_case = tdc_complete_case_analysis_robust(data, cox_covariates),
      augmented_complete_case =
        tdc_augmented_complete_case_analysis_robust(data, cox_covariates),
      nlac = tdc_nlac_analysis_robust(data, cox_covariates),
      iplw = tdc_ipw_analysis(data, cox_covariates, linkage_covariates)
    )
  } else {
    switch(
      method,
      oracle = ti_full_cohort_analysis(data, cox_covariates),
      complete_case = ti_complete_case_analysis(data, cox_covariates),
      augmented_complete_case =
        ti_augmented_complete_case_analysis(data, cox_covariates),
      nlac = ti_nlac_analysis(data, cox_covariates),
      iplw = ti_ipw_analysis(data, cox_covariates, linkage_covariates)
    )
  }

  if (identical(method, "iplw")) {
    coefficients <- result$beta_hat
    standard_errors <- result$beta_hat_se
    covariance <- result$beta_var / nrow(data$df)
    model <- NULL
  } else {
    coefficients <- stats::coef(result)
    covariance <- stats::vcov(result)
    standard_errors <- sqrt(diag(covariance))
    model <- result
  }

  structure(
    list(
      method = method,
      time_dependent = isTRUE(time_dependent),
      coefficients = coefficients,
      standard_errors = standard_errors,
      covariance = covariance,
      model = model,
      raw_result = result,
      call = match.call()
    ),
    class = "longterm_fit"
  )
}

#' Read the bundled synthetic PCPT dataset
#'
#' @param path Optional path to a compatible space-delimited dataset. By
#'   default, reads the synthetic dataset bundled with the package.
#'
#' @return A data frame with one row per synthetic participant.
#' @export
read_pcpt_synthetic <- function(path = NULL) {
  if (is.null(path)) {
    path <- system.file(
      "extdata", "PCPT_syntheticdata.txt", package = "longtermLinkage"
    )
  }
  if (!nzchar(path) || !file.exists(path)) {
    stop("The synthetic PCPT data file was not found.", call. = FALSE)
  }
  utils::read.table(path, header = TRUE, sep = " ", stringsAsFactors = TRUE)
}

#' Prepare PCPT data with time-dependent treatment effects
#'
#' @param data A subject-level data frame using the columns in
#'   `PCPT_syntheticdata.txt`.
#' @param change_points Two positive follow-up times defining treatment-effect
#'   change points. Defaults are 6.5 and 7.5 years expressed in days.
#'
#' @return A list with subject-level `df` and counting-process `tdata`.
#' @export
prepare_pcpt_linkage_data <- function(data, change_points = c(2373, 2738)) {
  required <- c(
    "id", "Y", "Delta", "Y_1", "Delta_1", "linkage", "txg",
    "censored_in_trial", "failed_in_trial"
  )
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop(
      "Missing required columns: ", paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.numeric(change_points) || length(change_points) != 2L ||
      any(!is.finite(change_points)) || any(change_points <= 0) ||
      change_points[1] >= change_points[2]) {
    stop("`change_points` must contain two increasing positive times.", call. = FALSE)
  }

  n <- nrow(data)
  first <- rep(change_points[1], n)
  second <- rep(change_points[2], n)
  tdata <- survival::tmerge(
    data, data, id = id, death = event(Y, Delta),
    T1 = tdc(first, txg)
  )
  tdata$T1[is.na(tdata$T1)] <- 0
  tdata <- survival::tmerge(
    tdata, data, id = id, T2 = tdc(second, txg)
  )
  tdata$T2[is.na(tdata$T2)] <- 0
  list(df = data, tdata = tdata)
}

#' Compare linkage methods on PCPT-format data
#'
#' @param data A list returned by `prepare_pcpt_linkage_data()`.
#' @param time_dependent Whether to include the two treatment-effect change
#'   terms `T1` and `T2`.
#' @param cox_covariates Optional Cox-model covariates.
#' @param linkage_covariates Optional linkage-model covariates.
#'
#' @return A list containing estimates, standard errors, confidence intervals,
#'   and fitted objects for IPLW, CC, CC+, and NLAC.
#' @export
compare_pcpt_methods <- function(
    data,
    time_dependent = TRUE,
    cox_covariates = NULL,
    linkage_covariates = c(
      "PSA", "racecode", "bmir", "famhx", "txg", "agebase",
      "EDUCATN", "MARISTAT", "EMPLSTAT", "PROFESN"
    )) {
  .validate_analysis_data(data, time_dependent)
  if (is.null(cox_covariates)) {
    cox_covariates <- if (isTRUE(time_dependent)) {
      c("PSA", "racecode", "bmir", "famhx", "txg", "agebase", "T1", "T2")
    } else {
      c("PSA", "racecode", "bmir", "famhx", "txg", "agebase")
    }
  }
  .check_covariates(cox_covariates, "cox_covariates")
  .check_covariates(linkage_covariates, "linkage_covariates")

  if (isTRUE(time_dependent)) {
    tdc_real_data_analysis(data, cv = cox_covariates, lv = linkage_covariates)
  } else {
    ti_real_data_analysis(data, cv = cox_covariates, lv = linkage_covariates)
  }
}

#' Run one simulation replicate from the paper
#'
#' @param n Positive integer sample size.
#' @param seed Positive integer random seed.
#' @param method One of `"oracle"`, `"complete_case"`,
#'   `"augmented_complete_case"`, `"nlac"`, `"iplw"`, or `"aiplw"`.
#' @param setting Integer from 1 through 7.
#' @param time_dependent Whether to use the time-dependent analysis. Defaults
#'   to `TRUE` for settings 1--5 and `FALSE` for settings 6--7.
#' @param n_boot Number of bootstrap samples for AIPLW.
#'
#' @return A list with coefficient estimates, standard errors, and coverage
#'   indicators for the known simulation truth.
#' @export
run_linkage_simulation <- function(
    n = 2000L,
    seed = 1L,
    method = c(
      "iplw", "nlac", "complete_case", "augmented_complete_case", "oracle",
      "aiplw"
    ),
    setting = 1L,
    time_dependent = setting <= 5L,
    n_boot = 100L) {
  method <- match.arg(method)
  .check_positive_integer(n, "n")
  .check_positive_integer(seed, "seed")
  .check_positive_integer(n_boot, "n_boot")
  setting <- as.integer(setting)
  if (length(setting) != 1L || is.na(setting) || !setting %in% 1:7) {
    stop("`setting` must be one integer from 1 through 7.", call. = FALSE)
  }
  if (setting %in% 6:7 && isTRUE(time_dependent)) {
    stop("Settings 6 and 7 support only time-independent analysis.", call. = FALSE)
  }
  if (setting == 4L && !isTRUE(time_dependent)) {
    stop("Setting 4 supports only time-dependent analysis.", call. = FALSE)
  }
  if (identical(method, "aiplw") && !setting %in% 6:7) {
    stop("AIPLW is implemented for settings 6 and 7.", call. = FALSE)
  }

  method_number <- match(
    method,
    c("oracle", "complete_case", "augmented_complete_case", "nlac", "iplw", "aiplw")
  )
  all_simulation(
    n = as.integer(n),
    seed = as.integer(seed),
    method_num = method_number,
    setting = setting,
    tdcAnalysis = isTRUE(time_dependent),
    robust = TRUE,
    n_BT = as.integer(n_boot)
  )
}

#' Methods for fitted long-term linkage models
#'
#' @param x,object A `longterm_fit` object.
#' @param parm Coefficients to include, by position or name.
#' @param level Confidence level.
#' @param ... Additional arguments, currently unused.
#'
#' @return `coef()` returns a named numeric vector, `vcov()` returns a matrix,
#'   and `confint()` returns a two-column matrix.
#' @name longterm_fit_methods
NULL

#' @rdname longterm_fit_methods
#' @export
print.longterm_fit <- function(x, ...) {
  cat("Long-term linkage Cox fit\n")
  cat("Method:", x$method, "\n")
  cat("Time-dependent:", x$time_dependent, "\n\n")
  table <- cbind(Estimate = x$coefficients, `Std. Error` = x$standard_errors)
  print(table)
  invisible(x)
}

#' @rdname longterm_fit_methods
#' @export
coef.longterm_fit <- function(object, ...) object$coefficients

#' @rdname longterm_fit_methods
#' @export
vcov.longterm_fit <- function(object, ...) object$covariance

#' @rdname longterm_fit_methods
#' @export
confint.longterm_fit <- function(object, parm, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be strictly between zero and one.", call. = FALSE)
  }
  estimate <- object$coefficients
  standard_error <- object$standard_errors
  if (!missing(parm)) {
    estimate <- estimate[parm]
    standard_error <- standard_error[parm]
  }
  critical <- stats::qnorm(1 - (1 - level) / 2)
  answer <- cbind(
    estimate - critical * standard_error,
    estimate + critical * standard_error
  )
  colnames(answer) <- paste0(format(c((1 - level) / 2, 1 - (1 - level) / 2) * 100), "%")
  answer
}

.check_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1 || x != as.integer(x)) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
}

.check_covariates <- function(x, name) {
  if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x))) {
    stop(sprintf("`%s` must be a non-empty character vector.", name), call. = FALSE)
  }
}

.validate_analysis_data <- function(data, time_dependent) {
  if (!is.list(data) || !is.data.frame(data$df)) {
    stop("`data` must be a list containing a subject-level `df` data frame.", call. = FALSE)
  }
  if (isTRUE(time_dependent) && !is.data.frame(data$tdata)) {
    stop("A time-dependent analysis requires a `tdata` data frame.", call. = FALSE)
  }
}
