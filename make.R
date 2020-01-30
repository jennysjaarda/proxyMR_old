#!/bin/bash
#SBATCH --partition sgg
#SBATCH --workdir /data/sgg2/jenny/projects/proxyMR
#SBATCH --job-name master
#SBATCH --output master.out
#SBATCH --mail-type=ALL                                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jennysjaarda@gmail.com                         # Where to send mail
#SBATCH --ntasks=1                                                 # Run on a single core


# make.R

# Master script for running proxyMR analyses.
# This project uses drake to manage workflows.
# For more information about drake, see
# https://ropensci.github.io/drake/

#################################
# Setup ----

# Set working directory
setwd(here::here())

# Load packages ----
source("code/packages.R")

# Load functions and plans  ----
source("code/settings.R")
source("code/functions.R")
source("code/plan.R")

if(file.exists("proxymr.log")){
  file.remove("proxymr.log")
}

# Run plans ----

partition = "sgg"

make(pre_pipeline,parallelism = "clustermq", console_log_file = "proxymr.log", cache_log_file = "cache_log.csv",
  memory_strategy = "lookahead", garbage_collection = TRUE, ## try with preclean ?
  jobs = 5, template = list(cpus = 1, partition = partition,
  log_file="/data/sgg2/jenny/projects/proxyMR/proxymr_%a_clustermq.out"))


make(full_proxymr,parallelism = "clustermq", console_log_file = "proxymr.log",  cache_log_file = "cache_log.csv",
  memory_strategy = "lookahead", garbage_collection = TRUE, ## try with preclean ?
  jobs = 94, template = list(cpus = 1, partition = partition,
  log_file="/data/sgg2/jenny/projects/proxyMR/proxymr_%a_clustermq.out"))

#
# make(run_mr, parallelism = "clustermq", console_log_file = "proxymr.log", cache_log_file = "cache_log.csv",
#   jobs = 75, template = list(cpus = 1, partition = partition,
#   log_file="/data/sgg2/jenny/projects/proxyMR/proxymr_%a_clustermq.out"))

#make(prep_shiny, console_log_file = "proxymr.log", cache_log_file = "cache_log.csv")


# read_csv("cache_log.csv", col_types = cols()) %>%
#   left_join(drake_cache_log(), by = "name") %>%
#   filter(hash.x != hash.y) %>%
#   dplyr::select(name, hash.x, hash.y, -type.x, -type.y)
