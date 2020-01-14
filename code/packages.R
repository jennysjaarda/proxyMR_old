library(processx)
library(pacman)
library(tidyverse)
library(drake)
library(ggplot2)
library(knitr)
library(purrr)
library(rlang)
library(tibble)
library(processx)
library(future)
library(future.batchtools)
library(pacman)
library(clustermq)
library(data.table)
library(devtools)
library(rmeta)
library(rslurm)
library(stringr)
library(bit64)
library(MASS)
library(qqman)
library(taRifx)
library(DataExplorer)
library(futile.logger)
library(tryCatchLog)
library(ggthemes)
library(readxl)
library(optparse)
library(reader)
library(TwoSampleMR)

install_github("ropensci/drake")
install.packages(c("DataExplorer","TwoSampleMR","bit64","clustermq","devtools","dplyr","futile.logger", "future", "future.batchtools","ggthemes","optparse","pacman","purr","qqman","reader","readxl","rmeta","rslurm","taRifx","tidyverse","tryCatchLog"))


The downloaded source packages are in
        ‘/tmp/Rtmp3psbDp/downloaded_packages’
Warning messages:
1: packages ‘TwoSampleMR’, ‘purr’ are not available (for R version 3.5.3)
2: In install.packages(c("DataExplorer", "TwoSampleMR", "bit64", "clustermq",  :
  installation of package ‘rmarkdown’ had non-zero exit status
3: In install.packages(c("DataExplorer", "TwoSampleMR", "bit64", "clustermq",  :
  installation of package ‘DataExplorer’ had non-zero exit status

devtools::install_github("boxuancui/DataExplorer")
