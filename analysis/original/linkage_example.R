# linkage_example.R
# This script runs simulation experiments and saves outputs to disk.

# Assumes the following file has already been sourced by run_all.R:
# source("linkage_functions.R")

# setwd first


# all simulation wrapper
simulation_all_wrapper <- function(start, n_rep, n, setting = 5, 
                                   tdcAnalysis = FALSE,
                                   method_num = 1,  
                                  output_dir = "./simulation_result/") {
  cat("Running simulation: method =", method_num, "| n =", n, "| setting =", setting, "\n")
  
  beta_1_arr = beta_2_arr = beta_3_arr = numeric(n_rep)
  cover_beta_1 = cover_beta_2 = cover_beta_3 = numeric(n_rep)
  beta_1_se = beta_2_se = beta_3_se = numeric(n_rep)
  
  seed_add_c = 1000
  start_index = (start - 1) * n_rep + 1
  end_index = start * n_rep
  i = start_index
  seed_add_times = 0
  
  while (i <= end_index) {
    if (i %% 20 == 0) cat("Iteration", i, "\n")
    
    skip_to_next = FALSE
    i_seed = i + seed_add_times * seed_add_c
    
    tryCatch({
      res <- all_simulation(n,  i_seed, method_num = method_num, 
                           setting = setting,
                           tdcAnalysis = tdcAnalysis,
                           robust = TRUE, n_BT = 100) 
    }, error = function(e) {
      cat("Error for seed", i_seed, ":", e$message, "\n")
      skip_to_next <<- TRUE
    })
    
    if (skip_to_next) {
      seed_add_times = seed_add_times + 1
      next
    }
    
    array_index = i - (start - 1) * n_rep
    beta_1_arr[array_index] = res$beta[1]
    beta_2_arr[array_index] = res$beta[2]
    beta_3_arr[array_index] = res$beta[3]
    cover_beta_1[array_index] = res$covered_beta_1
    cover_beta_2[array_index] = res$covered_beta_2
    cover_beta_3[array_index] = res$covered_beta_3
    beta_1_se[array_index] = res$beta_hat_se[1]
    beta_2_se[array_index] = res$beta_hat_se[2]
    beta_3_se[array_index] = res$beta_hat_se[3]
    i = i + 1
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  write(beta_1_arr, file = paste0(output_dir, "/beta_1_array_method_", method_num, "_start_", start, "_n_rep_", n_rep, "_setting_", setting, '.txt'))
  write(beta_2_arr, file = paste0(output_dir, "/beta_2_array_method_", method_num, "_start_", start, "_n_rep_", n_rep, "_setting_", setting, '.txt'))
  write(beta_3_arr, file = paste0(output_dir, "/beta_3_array_method_", method_num, "_start_", start, "_n_rep_", n_rep, "_setting_", setting, '.txt'))
  write(cover_beta_1, file = paste0(output_dir, "/cover_beta_1_method_", method_num, "_start_", start, "_n_rep_", n_rep, "_setting_", setting, '.txt'))
  write(cover_beta_2, file = paste0(output_dir, "/cover_beta_2_method_", method_num, "_start_", start, "_n_rep_", n_rep, "_setting_", setting, '.txt'))
  write(cover_beta_3, file = paste0(output_dir, "/cover_beta_3_method_", method_num, "_start_", start, "_n_rep_", n_rep, "_setting_", setting, '.txt'))
  write(beta_1_se, file = paste0(output_dir, "/beta_1_se_method_", method_num, "_start_", start, "_n_rep_", n_rep, "_setting_", setting, '.txt'))
  write(beta_2_se, file = paste0(output_dir, "/beta_2_se_method_", method_num, "_start_", start, "_n_rep_", n_rep, "_setting_", setting, '.txt'))
  write(beta_3_se, file = paste0(output_dir, "/beta_3_se_method_", method_num, "_start_", start, "_n_rep_", n_rep, "_setting_", setting, '.txt'))
  
  return(list(beta_1_arr = beta_1_arr, beta_2_arr = beta_2_arr, beta_3_arr = beta_3_arr,
              cover_beta_1 = cover_beta_1, cover_beta_2 = cover_beta_2, cover_beta_3 = cover_beta_3,
              beta_1_se = beta_1_se, beta_2_se = beta_2_se, beta_3_se = beta_3_se))
}
