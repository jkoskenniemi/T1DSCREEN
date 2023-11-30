
rm(list=ls())
#Get function import_data
source("code/01-Functions.R")
pacman::p_load(tidyverse, ieugwasr, TwoSampleMR)


#Names of available GWAS
names_gwas <- c("IL2RA", "IL2RB", "IL6R", "IL6ST", "JAK2", "JAK3", "TYK2")

#Save file paths
names_eqtl <- paste0("data/export_harmonization/", names_gwas, "_eqtl_TwoSampleMR.csv")
names_t1d <- paste0("data/export_harmonization/", names_gwas, "_T1D_TwoSampleMR.csv")

names_eqtl[4] <- "data/export_harmonization/IL6ST_prot_anno_TwoSampleMR.csv"
#Harmonize all GWAS of eqtl and t1d and save them to a list
harmonized_dfs <- map2(names_eqtl, names_t1d, import_data, has_beta = FALSE)

#Add names of dfs to list
names(harmonized_dfs)  <- names_gwas

#Extract LD matrices (this takes quite a while)
ld_matrices <- map(names_gwas, get_ld_matrix)

#Trim the LD matrices
trimmed_ld_matrices <- map2(ld_matrices, harmonized_dfs, trim_matrix)

#Rename LD matrices
names(trimmed_ld_matrices) <- paste0(names_gwas)

#Write LD matrices
map2(names_gwas, trimmed_ld_matrices, ~{
  assign(.x, .y, envir = globalenv())
  save(list = .x, file = paste0("data/export_LDmat/", .x, "_LDmat.RData"))
})

