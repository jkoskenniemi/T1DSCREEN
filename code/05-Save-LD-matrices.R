

rm(list = ls())
#Get function import_data
plinkref_path = "/home/jajoko/Documents/plinkref/EUR"

source("code/01-Functions.R")
pacman::p_load(tidyverse, ieugwasr, TwoSampleMR)


# 5-1 Get LD matrices around the targets----------------------------------------

# 5-1-1 eQTL data --------------------------------------------
#Names of available GWAS
names_gwas <- c("IL2RA", "IL2RB", "IL6R", "IL6ST", "JAK2", "JAK3", "TYK2")

#Save file paths
names_eqtl <- paste0("data/export_harmonization/",
                     names_gwas,
                     "_eqtl_TwoSampleMR.csv")
names_t1d <- paste0("data/export_harmonization/",
                    names_gwas,
                    "_T1D_TwoSampleMR.csv")

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
map2(names_gwas, trimmed_ld_matrices, ~ {
  assign(.x, .y, envir = globalenv())
  save(list = .x,
       file = paste0("data/export_LDmat/", paste0(.x, "_eqtl"), "_LDmat.rda"))
})

#5-1-2 pQTL data------------------------------

#Names of available GWAS
names_gwas <- c("IL6R", "IL6ST")

#Save file paths
names_pqtl <- paste0("data/export_harmonization/",
                     names_gwas,
                     "_prot_anno_TwoSampleMR.csv")
names_t1d <- paste0("data/export_harmonization/",
                    names_gwas,
                    "_T1D_TwoSampleMR.csv")

#Harmonize all GWAS of eqtl and t1d and save them to a list
harmonized_dfs <- map2(names_pqtl, names_t1d, import_data, has_beta = TRUE)

#Add names of dfs to list
names(harmonized_dfs)  <- names_gwas

#Extract LD matrices (this takes quite a while)
ld_matrices <- map(names_gwas, get_ld_matrix)

#Trim the LD matrices
trimmed_ld_matrices <- map2(ld_matrices, harmonized_dfs, trim_matrix)

#Rename LD matrices
names(trimmed_ld_matrices) <- paste0(names_gwas)

#Write LD matrices
map2(names_gwas, trimmed_ld_matrices, ~ {
  assign(.x, .y, envir = globalenv())
  save(list = .x,
       file = paste0("data/export_LDmat/", paste0(.x, "_pqtl"), "_LDmat.rda"))
})


#The following 2 sections (5-2 and 5-3) includes variants observed in analyses/IL2RA.Rmd
#and analyses/TYK2.Rmd so that those scripts can be run without running ieugwasr::ld_matrix()

#5-2 Get LD matrix within IL2RA that includes our top target and previously identified variants-----------


#Get LD matrix for our top variants (see IL2RA.Rmd and idependent variants
#identified in Robertson et al.variants near IL2RA (not included in the )

IL2RA_SNPs_LD <- ld_matrix(
  variants = c(
    "rs61839660",
    "rs35285258",
    "rs56179589",
    "rs6602437",
    "rs947474",
    "rs12722495"
  ),
  plink_bin = genetics.binaRies::get_plink_binary(),
  #this should work for most users, however for me genetics.binaRies() installs a faulty version of plink
  # plink_bin = "/bin/plink1.9",
  bfile = plinkref_path,
  #change this to your own plink location
  with_alleles = FALSE
)

save(IL2RA_SNPs_LD, file = "data/export_LDmat/IL2RA_tophits.rda")

#5-3 Get LD matrix within TYK2 that includes the previously selected interesting functional variants----------------

load("data/export_functional_variants/phenoscanner-entries.RDS")
functional_variants_TYK2 <- functional_variants$TYK2$results %>%
  select(rsid) %>%  distinct()

TYK2_SNPs_LD <- ld_matrix(
  variants = c(
    "rs34536443",
    #previously known top functional variant
    "rs12720356",
    "rs2304256",
    #our top functional variants, see analyses/TYK2.Rmd
    "rs144309607",
    "rs34725611"
  ),
  #top SNPs for TYK2 eQTL and T1D risk, see analyses/TYK2.Rmd
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = plinkref_path,
  with_alleles = FALSE
)

save(TYK2_SNPs_LD, file = "data/export_LDmat/TYK2_tophits.rda")
