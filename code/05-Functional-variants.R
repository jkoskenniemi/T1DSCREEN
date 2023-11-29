
#Load packages

#If run for the first time
# install.packages("devtools")
# library(devtools)
# install_github("phenoscanner/phenoscanner")
# library(phenoscanner)

#1 Obtain functional variants from Phenoscanner--------------------------------

# names_gwas <- c("IFNAR2", "IL2RA", "IL2RB",
#                 "IL12B", "IL2RG", "IL6R",
#                 "IL6ST", "JAK1", "JAK2",
#                 "JAK3", "TYK2")
#
# functional_variants <- map(names_gwas, ~phenoscanner(genequery = .x))
# names(functional_variants) <- names_gwas
#
# save(functional_variants, file = "data/export/phenoscanner-functional-variants.RDS")
load("data/import/phenoscanner-functional-variants.RDS")

#2 Turn results to a simple data frame-----------------------------------------

#Select only $results (i.e. flatten)
functional_variants_flattened <- sapply(functional_variants, function(x) x[2])

#Change to data.frame
df <- bind_rows(functional_variants_flattened)

#Select only diabetes variants
df %>%
  filter(str_detect(trait, "diabetes")) %>%
  select(gene, rsid, eur, consequence, protein_position, amino_acids, trait, efo, study, pmid)

filtered_missense <- df %>%
  filter(consequence == "missense") %>%
  filter(duplicated(rsid) == FALSE) %>%
  select(rsid, gene, hg38_coordinates, eur, consequence, protein_position, amino_acids)

# 3 Load QTL and T1D data on functional variants-------------------------------

#Names of available GWAS
names_gwas <- c("IFNAR2", "IL2RA", "IL2RB", "IL6R",
                "JAK1", "JAK2", "JAK3", "TYK2")
#IL6ST

names_eqtl <- paste0("data/export/", names_gwas, "_eqtl_TwoSampleMR.csv")
names_t1d <- paste0("data/export/", names_gwas, "_T1D_TwoSampleMR.csv")

#Add IL6ST
names_eqtl <- c(names_eqtl, "data/export/IL6ST_prot_anno_TwoSampleMR.csv")
names_t1d  <- c(names_t1d, "data/export/IL6ST_T1D_TwoSampleMR.csv")

harmonized_dfs <- map2(names_eqtl, names_t1d, import_data, has_beta = TRUE)

#flatten
funcvars_MR <- bind_rows(harmonized_dfs)
funcvars_MR <- filter(funcvars_MR, SNP %in% filtered_missense$rsid)

filtered_missense <- filtered_missense %>%
  select(rsid, gene, eur, protein_position, amino_acids, consequence) %>%
  rename(SNP = rsid)
funcvars_MR <- left_join(funcvars_MR, filtered_missense)

#4 Run MR-----------------------------------------------------------
funcvars_MR <- filter(funcvars_MR, pval.exposure < 1e-04)
MR_funcvars_MR <- mr_singlesnp(funcvars_MR)

MR_funcvars_MR <- MR_funcvars_MR %>%
  select(SNP, b, se, p) %>%
  filter(!is.na(b) & SNP != "All - Inverse variance weighted") %>%
  mutate(OR = exp(b),
         ci95_lo = exp(b) - 1.96*se,
         ci95_hi = exp(b) + 1.96*se)


#4 Table3 on functional variants--------------------------------------------

#Select only beta,ses and p-values of T1D
t1d <- map(names_t1d, import)
t1d <- bind_rows(t1d)

filtered_missense_t1d <- left_join(filtered_missense, t1d)
filtered_missense_t1d <- filtered_missense_t1d %>%
  select(gene, SNP, eaf, protein_position, amino_acids, hm_chrom, beta, se, pval) %>%
  rename(MAF = eaf, beta_t1d = beta, se_t1d = se, p_t1d = pval)

#Select betas, ses, and p-values of T1D and QTL-data
funcvars_MR <- funcvars_MR %>%
  select(SNP, gene, effect_allele.exposure, other_allele.exposure, eur,
         protein_position, amino_acids, beta.exposure, beta.outcome, pval.exposure, pval.outcome) %>%
  rename(EA = effect_allele.exposure, NEA = other_allele.exposure,
         AF = eur, beta_exp = beta.exposure, beta_t1d = beta.outcome,
         p_exp = pval.exposure, p_t1d = pval.outcome) %>%
  mutate(p_exp = signif(p_exp, 2), p_t1d = signif(p_t1d, 2)) %>%
  select(SNP, beta_exp, p_exp)

#select MR results
MR_funcvars_MR

#Join all three
filtered_missense_t1d %>%
  left_join(funcvars_MR, by = "SNP") %>%
  left_join(MR_funcvars_MR, by = "SNP") %>%
  select(gene, SNP, MAF, beta_t1d, se_t1d, p_t1d, beta_exp, p_exp, b, se, p, OR, ci95_lo, ci95_hi) %>%
  mutate(across(!c(gene, SNP), ~signif(.x, 2))) %>% #for nicer display
  write.csv("data/export/functional-variants-Table3.csv")

#expand
write.csv("data/export/functional_variants_t1d.csv")


