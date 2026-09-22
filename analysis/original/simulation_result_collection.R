loc = "./simulation_result/"

read_result <- function(setting, method_param, n_rep = 1000) {
  # Read simulation results
  #if (method_param == 6) {
  #  num_threads = 20
  #} else {
  #  num_threads = 1
  #}
  #n_rep = n_total / num_threads
  
  beta_1_arr = rep(0, n_rep)
  beta_2_arr = rep(0, n_rep)
  beta_3_arr = rep(0, n_rep)
  
  cover_beta_1 = rep(0, n_rep)
  cover_beta_2 = rep(0, n_rep)
  cover_beta_3 = rep(0, n_rep)
  
  beta_1_se = rep(0, n_rep)
  beta_2_se = rep(0, n_rep)
  beta_3_se = rep(0, n_rep)
  
  #for (i in 1 : num_threads) {
  i = 1
    beta_1_arr_filename = paste(
        loc, "beta_1_array_method_", method_param, "_start_", i, "_n_rep_", n_rep, 
        "_setting_", setting, '.txt', sep=""
      )
    beta_2_arr_filename = paste(
        loc, "beta_2_array_method_", method_param, "_start_", i, "_n_rep_", n_rep,
        "_setting_", setting, '.txt', sep=""
      )
    beta_3_arr_filename = paste(
        loc, "beta_3_array_method_", method_param, "_start_", i, "_n_rep_", n_rep, 
        "_setting_", setting, '.txt', sep=""
      )
    
    cover_beta_1_filename = paste(
      loc, "cover_beta_1_method_", method_param, "_start_", i, "_n_rep_", n_rep, 
      "_setting_", setting, '.txt', sep=""
    )
    cover_beta_2_filename = paste(
      loc, "cover_beta_2_method_", method_param, "_start_", i, "_n_rep_", n_rep, 
      "_setting_", setting, '.txt', sep=""
    )
    cover_beta_3_filename = paste(
      loc, "cover_beta_3_method_", method_param, "_start_", i, "_n_rep_", n_rep,
      "_setting_", setting, '.txt', sep=""
    )
    
    beta_1_se_filename = paste(
      loc, "beta_1_se_method_", method_param, "_start_", i, "_n_rep_", n_rep, 
      "_setting_", setting, '.txt', sep=""
    )
    beta_2_se_filename = paste(
      loc, "beta_2_se_method_", method_param, "_start_", i, "_n_rep_", n_rep, 
      "_setting_", setting, '.txt', sep=""
    )
    beta_3_se_filename = paste(
      loc, "beta_3_se_method_", method_param, "_start_", i, "_n_rep_", n_rep, 
      "_setting_", setting, '.txt', sep=""
    )
    
    start_index = (i - 1) * n_rep + 1
    end_index = i * n_rep
    beta_1_arr[start_index : end_index] = scan(beta_1_arr_filename)
    beta_2_arr[start_index : end_index] = scan(beta_2_arr_filename)
    beta_3_arr[start_index : end_index] = scan(beta_3_arr_filename)
    
    cover_beta_1[start_index : end_index] = scan(cover_beta_1_filename)
    cover_beta_2[start_index : end_index] = scan(cover_beta_2_filename)
    cover_beta_3[start_index : end_index] = scan(cover_beta_3_filename)
    
    beta_1_se[start_index : end_index] = scan(beta_1_se_filename)
    beta_2_se[start_index : end_index] = scan(beta_2_se_filename)
    beta_3_se[start_index : end_index] = scan(beta_3_se_filename)
  #}
  mean_beta_1 = mean(beta_1_arr)
  mean_beta_2 = mean(beta_2_arr)
  mean_beta_3 = mean(beta_3_arr)

  coverage_beta_1 = mean(cover_beta_1)
  coverage_beta_2 = mean(cover_beta_2)
  coverage_beta_3 = mean(cover_beta_3)

  mean_beta_1_se = mean(beta_1_se)
  mean_beta_2_se = mean(beta_2_se)
  mean_beta_3_se = mean(beta_3_se)

  beta_1_sd = sd(beta_1_arr)
  beta_2_sd = sd(beta_2_arr)
  beta_3_sd = sd(beta_3_arr)
  
  return(list("mean_beta_1"=mean_beta_1, "mean_beta_2"=mean_beta_2, "mean_beta_3"=mean_beta_3,
              "coverage_beta_1"=coverage_beta_1, "coverage_beta_2"=coverage_beta_2, "coverage_beta_3"=coverage_beta_3,
              "mean_beta_1_se"=mean_beta_1_se, "mean_beta_2_se"=mean_beta_2_se, "mean_beta_3_se"=mean_beta_3_se,
              "beta_1_sd"=beta_1_sd, "beta_2_sd"=beta_2_sd, "beta_3_sd"=beta_3_sd))
}

# Aggregates simulation results and prints summary metrics.

summarize_and_print_results <- function(setting, method_param, n_rep = 1000, 
                                        output_dir = "./simulation_result/") {
  source("simulation_result_collection.R")  # Load helper function
  
  result <- read_result(setting = setting, method_param = method_param, n_rep = n_rep)
  cat("Method", method_param, "| n =", n_rep, "\n")
  cat("Mean estimates: beta_1 =", result$mean_beta_1, 
      "| beta_2 =", result$mean_beta_2, 
      "| beta_3 =", result$mean_beta_3, "\n")
  cat("Coverage: beta_1 =", result$coverage_beta_1, 
      "| beta_2 =", result$coverage_beta_2, 
      "| beta_3 =", result$coverage_beta_3, "\n")
  cat("SE mean: beta_1 =", result$mean_beta_1_se, 
      "| beta_2 =", result$mean_beta_2_se, 
      "| beta_3 =", result$mean_beta_3_se, "\n")
  cat("Empirical SD: beta_1 =", result$beta_1_sd, 
      "| beta_2 =", result$beta_2_sd, 
      "| beta_3 =", result$beta_3_sd, "\n")
  return(result)
}
