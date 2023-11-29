#Build script

#if needed:
#install.packages("pacman")
#install.packages("devtools")
#devtools::install_github("phenoscanner")
#install_github("phenoscanner/phenoscanner")
pacman::p_load(data.table, dtplyr, dplyr, tidyverse, rio, here, TwoSampleMR, phenoscanner, phenoscanner)


source("code/01-Functions.R")
source("code/02-Gather.R")
source("code/03-Harmonize.R")
source("code/04-Save-LD-matrices.R")
source("code/05-Functional-variants.R")
