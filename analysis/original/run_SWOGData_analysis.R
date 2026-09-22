
#' Analysis of PCPT linked Data (using synthetic data as example) with Time-Dependent Covariates
#'
#' This script demonstrates how to analyze the PCPT linked data (including the synthetic data 
#' created from the PCPT linked dataset).
#' It fits survival models incorporating time-dependent treatment effects via change points
#' using methods such as IPW, CC, ACC, and NLAC.
#'
#' ## Steps:
#' 1. Read synthetic data from file
#' 2. Construct time-dependent covariates (`T1`, `T2`) using `tmerge()`
#' 3. Prepare covariates for Cox and linkage models
#' 4. Run the main analysis using `tdc_real_data_analysis()` and `ti_real_data_analysis()`
#'

source('SWOGData_analysis_functions.R')
# ---- Load Synthetic PCPT Data ----
## Read the synthetic dataset with survival and covariate information
syntheticData = read.table(file = './_SyntheticDataSet/PCPT_syntheticdata.txt', header = TRUE, sep = ' ')

# ---- Define Change Points and Create Time-Dependent Covariates ----

## First change point: 6.5 years (~2373 days)
first_change_point = 2373 # first change point at 6.5 years.
t_1 = rep(first_change_point, nrow(syntheticData))
tdata = tmerge(syntheticData, syntheticData, id=id, death = event(Y, Delta),
               T1 = tdc(t_1, txg))
tdata$T1[is.na(tdata$T1)] = 0

## Second change point: 7.5 years (~2738 days)
second_change_point = 2738 # second change point at 7.5 years. 
t_2 = rep(second_change_point, nrow(syntheticData))
tdata = tmerge(tdata, syntheticData, id=id, T2 = tdc(t_2, txg))
tdata$T2[is.na(tdata$T2)] = 0 


# ---- Wrap Data into Required Format ----
syntheticanalysisdata = list("tdata"=tdata, "df"=syntheticData)


# ---- Specify Covariates for Models ----
## Time-invariant covariates for Cox model
cox_covariates <- c("PSA", "racecode", "bmir", "famhx",
                    "txg", "agebase")

## Additional covariates for linkage (logistic regression) model
additional_logistic_covariates <- c("EDUCATN", "MARISTAT", "EMPLSTAT", "PROFESN")

## Cox covariates including time-dependent components
td_cox_covariates = c(cox_covariates, "T1", "T2") # add time-dependent covariates


# ---- Run Analysis with Robust SEs (time independent and dependent covariates)----
# produce the results from sysnthetic data similar to the results in Table 3 in main paper
tdc_result = tdc_real_data_analysis(data = syntheticanalysisdata, 
                                   cv=td_cox_covariates,
                                   lv=c(cox_covariates, additional_logistic_covariates))

# ---- Run Analysis time independent covariates)---- 
# produce the results from sysnthetic data similar to the results in Table 7 in supplementary material
ti_result = ti_real_data_analysis(data = syntheticanalysisdata, 
                                    cv=cox_covariates,
                                    lv=c(cox_covariates, additional_logistic_covariates))
