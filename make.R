# make.R

# Master script for running PSYMETAB analyses.
# This project uses drake to manage workflows.
# For more information about drake, see
# https://ropensci.github.io/drake/

# Setup ----

# Set working directory
setwd(here::here())

# Load packages ----
source("code/packages.R")

# Load functions and plans  ----
source("code/settings.R")
source("code/functions.R")
source("code/plan.R")



make(pre_pipeline,parallelism = "clustermq",console_log_file = "proxymr.log", cache_log_file = "cache_log.csv",
  memory_strategy = "lookahead", garbage_collection = TRUE, ## try with preclean ?
  jobs = 5, template = list(cpus = 1, partition = "sgg",
  log_file="/data/sgg2/jenny/projects/proxyMR/proxymr_%a_clustermq.out"))


make(full_proxymr,parallelism = "clustermq",console_log_file = "proxymr.log",  cache_log_file = "cache_log.csv",
  memory_strategy = "lookahead", garbage_collection = TRUE, ## try with preclean ?
  jobs = 120, template = list(cpus = 1, partition = "sgg",
  log_file="/data/sgg2/jenny/projects/proxyMR/proxymr_%a_clustermq.out"))
