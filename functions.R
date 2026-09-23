###############################################################################
### Function ##################################################################
### - Section 1: Function for analyse
###############################################################################

# load packages
library(tidyverse)

library(tidyverse)

# extract result out of cox model
extract_model_results <- function(model, model_name, true_beta) {
  
  # model summary
  mod_sum <- summary(model)
  coef_mat <- mod_sum$coefficients
  
  # extract values
  log_hr <- coef_mat["arm", "coef"]
  p_val <- coef_mat["arm", "Pr(>|z|)"]
  
  # check for robust se
  se_col <- ifelse("robust se" %in% colnames(coef_mat), "robust se", "se(coef)")
  robust_se <- coef_mat["arm", se_col]
  
  # summary in tibble
  tibble(
    Model = model_name,
    True_Beta = true_beta,
    Log_HR = round(log_hr, 4),
    Hazard_Ratio = round(exp(log_hr), 2),  
    Robust_SE = round(robust_se, 4),
    P_Value = round(p_val, 4),
    Significant = ifelse(p_val < 0.05, "Yes", "No"),
    Bias = round(log_hr - true_beta, 4)
  )
}