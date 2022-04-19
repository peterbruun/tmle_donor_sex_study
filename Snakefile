# Snakemake run file for analysis


import numpy as np
import pandas as pd
import datetime
pd.set_option('max_info_columns', 10**10)
pd.set_option('display.max_columns', 1000)
pd.set_option('max_info_columns', 10**10)
pd.set_option('max_info_rows', 10**10)

workdir: "/users/data/projects/deep_phenotyping/TMLE_donor_sex_and_age/"


# Recipient stratification
RECIPIENT_SEX = ["M","F"]

# Sex analyses
TREATMENT_SEX = ["male_donors","female_donors","natural_course"]


# Follow_up
FOLLOW_UP = ["28"]

# Fatnode memory request
memory = 1024*120

# Threads for multiprocessing
NTHREADS = 40



###### Main rule #######

rule all:
	input:		
		# Donor sex analyses
		# 5-fold cross-validation
		expand("results/donor_sex/2022-02-20/{recipient_sex}_recipients/mortality/{treatment}_{follow_up}days_kInf_trt05_5fold.rds", 
			recipient_sex = RECIPIENT_SEX, treatment = TREATMENT_SEX, follow_up = FOLLOW_UP),
		"results/donor_sex/2022-02-20/analysis_code.R",


		# Input files
		# 28-days follow-up available
		longitudinal_data_mortality_28 = "/data/projects/deep_phenotyping/transfusions_bth_simple/data/processed/data_longitudinal_inclAcute_mortality_28days.tsv",


## ------------------------------------------ RULES ------------------------------------------ ##




## ----------------- DONOR SEX ------------------------ ##



rule model2_save_code_donor_sex:
	input:
		code = "scripts/donor_sex/sm_hal_lmtp_donor_sex.R",
	output:
		code = "results/donor_sex/2022-02-20/analysis_code.R",
	resources:
		tmin = 10,
		mem_mb = 1024*2,
	priority: 2
	shell:
		"""
		cp {input.code} {output.code}
		"""

rule model2_lmtp_donor_sex:
	input:
		data = rules.all.input.longitudinal_data_mortality_28,
	output:
		survival = "results/donor_sex/2022-02-20/{recipient_sex}_recipients/mortality/{treatment}_{follow_up}days_kInf_trt05_5fold.rds"
	params:
		exposure = "male_donor",
		follow_up = "{follow_up}",
		recipient_sex = "{recipient_sex}",
		intervention = "{treatment}",
	resources:
		tmin = 60*24*7,
		mem_mb = memory,
	threads: NTHREADS,
	priority: 1,
	run:
		# Run analysis
		lmtp_model = [
			'Rscript scripts/donor_sex/sm_lmtp_donor_sex.R',
			'{input.data}',
			'{params.exposure}',
			'{params.follow_up}',
			'{params.recipient_sex}',
			'{output.survival}',
			"{params.intervention}",
			"{threads}"]
		shell(' '.join(lmtp_model))







