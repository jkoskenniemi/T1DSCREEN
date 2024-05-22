#Build script

#You need Plink installed in your system to run these analyses
#Set path to Plink here:
plink_path = "~/Documents/plinkref/EUR"

#if any of these tools is not installed:
#install.packages("pacman")
install.packages("devtools")
devtools::install_github("phenoscanner/phenoscanner")
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("MRCIEU/genetics.binaRies")

#Load packages
pacman::p_load(data.table, dtplyr, dplyr, tidyverse, rio, here, TwoSampleMR,
               phenoscanner, phenoscanner, purrr, ieugwasr, genetics.binaRies, coloc)



source("code/01-Functions.R")
source("code/02-Gather.R")
source("code/03-Harmonize.R")
source("code/04-Functional-variants.R")
source("code/05-Save-LD-matrices.R")
