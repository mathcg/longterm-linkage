library(longtermLinkage)

pcpt <- read_pcpt_synthetic()
data <- prepare_pcpt_linkage_data(pcpt)
results <- compare_pcpt_methods(data, time_dependent = TRUE)

results$beta_hat
results$beta_hat_se
