#####################################################################################
## Author: Jenny Sjaarda
## Project: MR_Shared_Environment
##
## Time-stamp: <[settings.r] by JS 2019-04-02 10:42:03 CEST>
##
## Description:
##
##
## History:
##
#####################################################################################

### Set data and project directories

project_dir <- "/data/sgg2/jenny/projects/MR_Shared_Environment/"
SGG_generic <- "/data/sgg2/jenny/SGG_generic/"

source(paste0(SGG_generic,"/scripts/settings.r"))
source(paste0(project_dir,"/scripts/functions.r"))

traits <- c()
MR_method_list <- c("mr_wald_ratio","mr_ivw","mr_ivw_fe","mr_egger_regression","mr_weighted_median")

### define variables
IV_threshold <- 5e-08
prune_threshold <- 0.001
GRS_thresholds <- c(0.1,0.01,0.001)
num_couple_length_bins <- 10

num_IVs_threshold <- 5
household_correlation_threshold <-0.1
irnt=TRUE

phesant_directory <- read.table(paste0(UKBB_processed,"/PHESANT/","PHESANT_file_directory.txt"), header=T)

# Data-Field 6141
# Description:	How are people in household related to participant
# limit to only individuals who responded YES to husband, wife, or partner

relatedness_field <- "6141_1"
relations_file <- as.character(phesant_directory[which(phesant_directory[,2]==relatedness_field),"File"])

# name of first phesant file
first_phesant_file <- as.character(phesant_directory[1,"File"])



time_at_address_file <- "/data/sgg2/jenny/data/UKBB_processed/PHESANT/ukb31459/bin1/out_bin1..tsv"
time_at_address_raw_file <- "/data/sgg2/jenny/data/UKBB_raw/pheno/ukb31459.csv"
time_at_address_field <- "699"
time_at_address_raw_field <- "699-0.0"
## functions to fit:
