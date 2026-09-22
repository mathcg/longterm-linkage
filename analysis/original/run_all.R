# run_all.R
# Master script to run the full simulation pipeline.
# Steps:
# 1. Load function definitions
# 2. Run simulation experiments (may take time)
# 3. Aggregate and summarize results

## this set of simulation will reproduce Table 1 results in the main paper
## (Simulation results when Cox model is correctly specified when n = 2000.)

## method_param = 1, Full cohort analysis (oracle method)
## method_param = 2, Complete case analysis
## method_param = 3, augmented complete case analysis
## method_param = 4, NLAC analysis
## method_param = 5, IPLW analysis 
## method_param = 6, AIPLW analysis 

## setting = 1, linkage is LCAR, time independent & dependent covariates
## setting = 2, linkage is CLAR, time independent & dependent covariates
## setting = 3, linkage is LNAR(\tilde T), time independent & dependent covariates
## setting = 4, linkage is LNAR(C_2), time independent & dependent covariates
## setting = 6, linkage is MAR, time independent & dependent covariates
## setting = 6, linkage is MAR, time independent covariates only
## setting = 7, linkage is CLAR, time independent covariates only

## tdcAnalysis = FALSE, run time independent analysis methods
## tdcAnalysis = TRUE, run time independent analysis methods

## To reproduce Table 1 in main paper, set setting = 1 to 4, method_param = 1 to 5, tdcAnalysis = TRUE
## To reproduce Table 2 in main paper, set setting = 1 to 3, method_param = 1 to 5, tdcAnalysis = FALSE
## To reproduce Table 1 in suppl, set setting = 5, method_param = 1 to 5, tdcAnalysis = TRUE
## To reproduce Table 2 in suppl, set setting = 5, method_param = 1 to 5, tdcAnalysis = FALSE
## To reproduce Table 3 in suppl, set setting = 6 and 7, method_param = 1 to 6, tdcAnalysis = FALSE

## setwd first

cat("Step 1: Loading core functions...\n")
source("linkage_functions.R")

cat("Step 2: Running simulations...\n")
# Set your parameters here
start <- 1
n_rep <- 1000 
setting <- 1
n <- 2000
method_param <- 4 # Change this to desired method (1 to 5)
tdcAnalysis <- TRUE
source("linkage_example.R") 
all_result <- simulation_all_wrapper(start, n_rep, n, setting = setting, 
                                     tdcAnalysis = tdcAnalysis,
                                   method_num = method_param, 
                                   output_dir = "./simulation_result/") 

cat("✅ Simulation complete. Intermediate results saved.\n\n")

cat("Step 3: Summarizing results...\n")
source("simulation_result_collection.R")
summarize_and_print_results(setting, method_param, n_rep = n_rep, 
                            output_dir = "./simulation_result/")
cat("\n🎯 Summary complete. See printed output for means, SEs, and coverage.\n")
