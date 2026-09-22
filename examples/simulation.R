library(longtermLinkage)

data <- simulate_linkage_data(n = 1000, setting = 1, seed = 2026)

fit <- fit_linkage_cox(
  data,
  method = "iplw",
  cox_covariates = c("Z_1", "Z_2", "Z_3"),
  linkage_covariates = c("Z_1", "Z_2"),
  time_dependent = TRUE
)

print(fit)
confint(fit)
