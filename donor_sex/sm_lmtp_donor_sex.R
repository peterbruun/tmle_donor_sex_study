#############################
# TMLE for donor sex

args <- commandArgs(trailingOnly=TRUE)

library("future")
library(SuperLearner)
library(lmtp)
library(hal9001)
library(tidyverse)
library(progressr)
library(broom)

setwd("/users/data/projects/deep_phenotyping/TMLE_donor_sex_and_age/")

# Testing
# infile = "/data/projects/deep_phenotyping/transfusions_bth_simple/data/processed/data_longitudinal_inclAcute_mortality.tsv"
# expo_var = "male_donor"
# follow_up = 28
# recipient_sex = "M"
# threads = 39
# outfile = "test.rds"



# Running
infile = args[1]
expo_var = as.character(args[2])
follow_up = as.numeric(args[3])
recipient_sex = as.character(args[4])
outfile = as.character(args[5])
intervention = as.character(args[6])
threads = as.numeric(args[7])



# 28 days
data_bootstrap = read_tsv(infile)


# Limit time wise
data_bootstrap <- data_bootstrap %>% filter(time < follow_up + 1)

# Limit to recipient sex
data_bootstrap <- data_bootstrap %>% filter(patient_cpr_sex == recipient_sex)

# arrange patients according to time
data_bootstrap <- data_bootstrap %>% arrange(patient_cpr_enc,time)

# Define exposure variable
data_bootstrap["trt"] <- data_bootstrap[expo_var]

# count variable
#data_bootstrap["exposure"] <- data_bootstrap["male_donor"]
# Set exposure to ratio between male/n_trans and -1 on days where no transfusions are received.
#data_bootstrap <- data_bootstrap %>% mutate(exposure = if_else(n_transfusions>0, male_donor/n_transfusions, -1))
# Set to 50/50 on days where no transfusions where received
data_bootstrap <- data_bootstrap %>% mutate(exposure = if_else(n_transfusions>0, trt/n_transfusions, 0.5))


# Make total_transfusions not included transfusions on this day k
data_bootstrap["lag1_total_transfusions"] = data_bootstrap["total_transfusions"] - data_bootstrap["n_transfusions"]
data_bootstrap <- data_bootstrap %>% select(-c(total_transfusions))


# data_bootstrap <- data_bootstrap %>% 
#                    dplyr::select(c(patient_cpr_enc,patient_cpr_sex,event,time,exposure,male_donor,n_transfusions,lag1_total_transfusions,))

# Select columns to use
data_bootstrap <- data_bootstrap %>% 
                   dplyr::select(c(patient_cpr_enc,event,time,exposure,donor_age_below_40,n_transfusions,lag1_total_transfusions,
                            patient_charlson_score,patient_age_trans,Trans_year,sin_month,cos_month,
                            AB0_patient,Rhesus_patient,hospital_navn))



#data_bootstrap <- model.matrix(data_bootstrap)

#one approach it to index with the $ sign and the as.factor function
#AB0_patient,Rhesus_patient,hospital_navn
data_bootstrap$AB0_patient <- as.factor(data_bootstrap$AB0_patient)
data_bootstrap$Rhesus_patient <- as.factor(data_bootstrap$Rhesus_patient)
data_bootstrap$hospital_navn <- as.factor(data_bootstrap$hospital_navn)
data_bootstrap$event <- as.integer(data_bootstrap$event)

# convert event to integer
data_bootstrap$event <- as.integer(data_bootstrap$event)

# Add censoring variable
data_bootstrap["C"] <- 1
# wide format
data_bootstrap <- data_bootstrap %>% 
										pivot_wider(names_from = time, 
																values_from = c(event,exposure,donor_age_below_40,
																	n_transfusions,lag1_total_transfusions,
                            			patient_charlson_score,Trans_year,sin_month,cos_month,C,
                            			hospital_navn))

# Rename patient id to rowname for memory efficiency
#data_bootstrap$patient_cpr_enc <- rownames(data_bootstrap)
data_bootstrap <- data_bootstrap %>% select(-c(patient_cpr_enc))

# If year and month should not be time-varying
data_bootstrap <- data_bootstrap %>% select(-c(
														 paste0("Trans_year_", 1:follow_up),
														 paste0("sin_month_", 1:follow_up),
														 paste0("cos_month_", 1:follow_up)))


Y <- paste0("event_", 0:follow_up)
data_bootstrap <- event_locf(data_bootstrap, Y)

#data_bootstrap %>% filter(event_3 == 1) %>% select(Y)

# Carry event forward - alternative version
# Y <- paste0("event_", 0:follow_up)
# data_bootstrap <- data_bootstrap %>%
#   mutate_at(Y, ~replace_na(., 1))


## --- TMLE --- ##

# Set treatment to be equal to L at time t
male_donors_only <- function(data, trt) {
  return(ifelse(
    (data[[sub("exposure", "n_transfusions", trt)]] > 0) & 
    !(is.nan(data[[sub("exposure", "n_transfusions", trt)]])), # If transfusions are received
    1, # set all to male
    (data[[trt]]) # otherwise do nothing
  ) )
}

# Set treatment to be equal to L at time t
female_donors_only <- function(data, trt) {
  return(ifelse(
    (data[[sub("exposure", "n_transfusions", trt)]] > 0) & 
    !(is.nan(data[[sub("exposure", "n_transfusions", trt)]])), # If transfusions are received
    0, # set all to female
    (data[[trt]]) # otherwise do nothing
  ) )
}

# Test treatment strategies
# test <- male_donors_only(data_male_recipients,"exposure_1")
# test2 <- data_male_recipients %>% select(c(exposure_1,n_transfusions_1))
# test2["updated_treat"] <- test


# For count treatment
# male_donors_only <- function(data, trt) {
#   (data[[trt]] = data[[sub("exposure", "n_transfusions", trt)]])
# }
# female_donors_only <- function(data, trt) {
#   (data[[trt]] = 0)
# }

# Define covariates
# Time-varying covariates
L = list() # Make an empty list to save output in
for (i in 1:follow_up+1) {

		# Non time-varying year and month
		columns <- c("donor_age_below_40_XX","n_transfusions_XX","lag1_total_transfusions_XX","patient_charlson_score_XX",
			"hospital_navn_XX")
		# Time-varying year and month    
		#columns <- c("donor_age_below_40_XX","n_transfusions_XX","lag1_total_transfusions_XX","patient_charlson_score_XX","Trans_year_XX","sin_month_XX","cos_month_XX","hospital_navn_XX")
    L[[i]] = str_replace(columns,"XX",as.character(i-1))
}

# Treatment
A <- paste0("exposure_", 0:follow_up)
# Outcome
Y <- paste0("event_", 0:follow_up)
# Censoring
C <- paste0("C_", 0:follow_up)
# Baseline covariates
W <- c("AB0_patient","Rhesus_patient","patient_age_trans",
			 "Trans_year_0","sin_month_0","cos_month_0")
# Time-varying year and month 
#W <- c("AB0_patient","Rhesus_patient","patient_age_trans")
tune = list(ntrees = c(1000), max_depth = c(2, 4),
            shrinkage = c(0.1, 0.01),minobspernode = c(10),nthread = c(threads-3))
xgb_grid = create.SL.xgboost(tune = tune)
elastic_net = create.Learner("SL.glmnet", tune = list(alpha = c(0,0.5,1)))

# glm.interaction takes way to long!
#lrnrs <- c("SL.glm","SL.glm.interaction","SL.glmnet",elastic_net$names,"SL.earth",rf_grid$names, xgb_grid$names)


lrnrs <- c("SL.glm","SL.glmnet","SL.earth",xgb_grid$names)



## --- TMLE --- ##
handlers(global = TRUE)


# Multi-proccessing
# Multi-proccessing
options(mc.cores = threads)
#plan(multisession,workers = threads)
plan(multicore,workers = threads)
set.seed(42, "L'Ecuyer-CMRG")

if (intervention == "male_donors"){

results <- lmtp_tmle(
			data = data_bootstrap, 
			trt = A, 
			outcome = Y, 
			baseline = W, 
			time_vary = L, 
			cens = C,
			k = Inf,
			intervention_type = "dynamic", 
			shift = male_donors_only, 
			outcome_type = "survival",
			learners_trt = lrnrs, 
      learners_outcome = lrnrs,
			folds = 5, 
			.SL_folds = 5)

} else if (intervention == "female_donors"){

results <- lmtp_tmle(
			data = data_bootstrap, 
			trt = A, 
			outcome = Y, 
			baseline = W, 
			time_vary = L, 
			cens = C,
			k = Inf, 
			intervention_type = "dynamic", 
			shift = female_donors_only, 
			outcome_type = "survival",
			learners_trt = lrnrs, 
      learners_outcome = lrnrs,
			folds = 5, 
			.SL_folds = 5)

} else if (intervention == "natural_course"){

results <- lmtp_tmle(
			data = data_bootstrap, 
			trt = A, 
			outcome = Y, 
			baseline = W, 
			time_vary = L, 
			cens = C,
			k = Inf,
			intervention_type = "dynamic", 
			shift = NULL, 
			outcome_type = "survival",
			learners_trt = lrnrs, 
      learners_outcome = lrnrs,
			folds = 5, 
			.SL_folds = 5)

}

# Estimate survival probability
# Save results
saveRDS(results, outfile)







