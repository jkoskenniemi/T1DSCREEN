#Load package

source("code/01-Functions.R")

# This script merges the annotation and pQTL/eQTL data.
#
# GWAS sumstat of type 1 diabetes was obtained from the study by Chiou et al.
# (Nature 2021, https://doi.org/10.1038/s41586-021-03552-w). GWAS sumstat for
# protein levels of IL2RG was obtained from the study by Ferkingstad et al.
# (Nature 2021, https://doi.org/10.1038/s41588-021-00978-w).
#
# Variants within 1 megabase before the start site of IL2RG and 1 megabase
# after the stop site are included. See '/code/gather.R' on how the the
# variants within this area were filtered in each GWAS. The gather-script
# was run in the high performance computing cluster using the
# 'run_gather.sh' script. The output of '/code/Gather.R' is saved in
# '/data/cis_sumstats/', and for each possible target gene
# there are three files:
#
# 1. *_T1D_data.rds: GWAS sumstat for the risk of T1D (Chiou et al. 2021)
# 2. *_prot_data.rds: GWAS sumstat of serum protein level (Ferkingstad et al. 2021)
# 3. *_anno_data.rds: further annotation for GWAS of serum protein level
# Where * is a gene of interest (IL2RA, IL2RG, IL2RG, IL6ST, IL6R, etc.)
#
# except for IL6R for which
# IL6R_prot_data.rds is substituted with IL6R_crp_data.rds

# 1 pQTL data------------------------------------------------------------------

## 1.1 Import pQTL data----------------------------------------------------------

#Load GWAS sumstat for IL6ST
IL6ST_prot <- import("data/export_cis_sumstats/IL6ST_prot.rds")
IL6ST_anno <- import("data/export_cis_sumstats/IL6ST_anno.rds")

#Load these for IL2RA, IL6R and TYK2 (again, reviewer's idea)
IL6R_prot <- import("data/export_cis_sumstats/IL6R_prot.rds")
IL6R_anno <- import("data/export_cis_sumstats/IL6R_anno.rds")
IL2RA_prot <- import("data/export_cis_sumstats/IL2RA_prot.rds")
IL2RA_anno <- import("data/export_cis_sumstats/IL2RA_anno.rds")
TYK2_prot <- import("data/export_cis_sumstats/TYK2_prot.rds")
TYK2_anno <- import("data/export_cis_sumstats/TYK2_anno.rds")

## 1.2 Merge protein data with annotations--------------------------------------

# Ferkingstad et al. state in their README file
#
#### A subset of the variants in the summary statistics files should be excluded
#### due to quality issues. These variants are listed in a separate Excluded variants
#### file (assocvariants.excluded.txt.gz). The Extra annotation file
#### (assocvariants.annotated.txt.gz) does not include these variants.
#
# Therefore in the join process, only the data available in the annotation (IL2RG_anno)
# is used (hence, left_join)

IL6ST_prot_anno <- left_join(IL6ST_prot, IL6ST_anno, by = c("Name", "Chrom", "Pos"))
IL6R_prot_anno <- left_join(IL6R_prot, IL6R_anno, by = c("Name", "Chrom", "Pos"))
IL2RA_prot_anno <- left_join(IL2RA_prot, IL2RA_anno, by = c("Name", "Chrom", "Pos"))
TYK2_prot_anno <- left_join(TYK2_prot, TYK2_anno, by = c("Name", "Chrom", "Pos"))
rm(IL6ST_anno, IL6ST_prot, IL6R_anno, IL6R_prot, IL2RA_anno, IL2RA_prot,
   TYK2_anno, TYK2_prot)

## 1.3 Sort conflicts-----------------------------------------------------------

#Further in the README file of Ferkingstad et al.
#
#### Also, the summary statistics files sometimes incorrectly have effectAllele=otherAllele for multiallelic variants.
#### In these cases the effectAllele is correct, but the otherAllele should be '!', meaning that the effectAllele is tested against the other (two or more) alleles (using the '!' sign as shorthand for 'not effectAllele').
#### This has been corrected in the file Extra annotation file (assocvariants.annotated.txt.gz).
#
#Therefore, multiallelic variants were removed.

filter_multiallelic <- function(x) filter(x, otherAllele.y != "!")
IL6ST_prot_anno <- filter_multiallelic(IL6ST_prot_anno )
TYK2_prot_anno <- filter_multiallelic(TYK2_prot_anno )
IL2RA_prot_anno <- filter_multiallelic(IL2RA_prot_anno )
IL6R_prot_anno <- filter_multiallelic(IL6R_prot_anno )

# After removal of multiallelic variants there were no conflicts
# between effectAllele.x and effectAllele.y or between
# otherAllele.x and otherAllele.y

#List those with conflicts in effectAllele or otherAllele
check_conflicts_ea_oa(IL6ST_prot_anno) #none
check_conflicts_ea_oa(IL6R_prot_anno) #none
check_conflicts_ea_oa(IL2RA_prot_anno) #none
check_conflicts_ea_oa(TYK2_prot_anno) #none

# After this, the remaining discrepancies between the variables
# obtained from annotation and protein data files were because those
# variants without an rsid were listed as "." in the annotation file
# and as "NA" in the main GWAS sumstat file

#List subjects which had conflicting rsid information and show
#which entries there are for these


#All the conflicts in rsid.x and rsid.y are because rsid.x is NA and rsid.y is "."
check_conflicts_rsid(IL6ST_prot_anno)
check_conflicts_rsid(TYK2_prot_anno)
check_conflicts_rsid(IL2RA_prot_anno)
check_conflicts_rsid(IL6R_prot_anno)

# Thus, only those values variables with .x in their name are preserved.
# Here, we also remove the remaining multiallelic SNPs
# (where effect allele equals the other allele.)


#Discard all variables that end with an ".y"
IL6ST_prot_anno <- IL6ST_prot_anno %>% select(!ends_with(".y"))
IL6R_prot_anno <- IL6R_prot_anno %>% select(!ends_with(".y"))
IL2RA_prot_anno <- IL2RA_prot_anno %>% select(!ends_with(".y"))
TYK2_prot_anno <- TYK2_prot_anno %>% select(!ends_with(".y"))

#Remove ".x" from variable names
IL6ST_prot_anno <- IL6ST_prot_anno %>% rename_with(~gsub(".x", "", .x, fixed = TRUE))
TYK2_prot_anno <- TYK2_prot_anno %>% rename_with(~gsub(".x", "", .x, fixed = TRUE))
IL2RA_prot_anno <- IL2RA_prot_anno %>% rename_with(~gsub(".x", "", .x, fixed = TRUE))
IL6R_prot_anno <- IL6R_prot_anno %>% rename_with(~gsub(".x", "", .x, fixed = TRUE))

#Sanity checks
colnames(IL6ST_prot_anno) #all colnames are clean (n colname.x or colname.y)
colnames(IL6R_prot_anno) #all colnames are clean (n colname.x or colname.y)
colnames(IL2RA_prot_anno) #all colnames are clean (n colname.x or colname.y)
colnames(TYK2_prot_anno) #all colnames are clean (n colname.x or colname.y)

# Entries without rsids are removed
IL6ST_prot_anno <- filter(IL6ST_prot_anno, !is.na(rsids))
TYK2_prot_anno <- filter(TYK2_prot_anno, !is.na(rsids))
IL6R_prot_anno <- filter(IL6R_prot_anno, !is.na(rsids))
IL2RA_prot_anno <- filter(IL2RA_prot_anno, !is.na(rsids))

## 1.4 Write data in TwoSampleMR format-----------------------------------------

# Data are written in TwoSampleMR format so that they can be easily opened
# in the main analyses scripts.

#Rename
IL6ST_prot_anno <- rename_for_TwoSampleMR_e(IL6ST_prot_anno)
IL2RA_prot_anno <- rename_for_TwoSampleMR_e(IL2RA_prot_anno)
IL6R_prot_anno <- rename_for_TwoSampleMR_e(IL6R_prot_anno)
TYK2_prot_anno <- rename_for_TwoSampleMR_e(TYK2_prot_anno)

# pqtl_names <- c("IL6ST", "IL6R", "IL2RA", "TYK2")
#Add phenotype information
# pqtl_phenotype_names <- paste0("Serum ", pqtl_names, " level")
IL6ST_prot_anno <- IL6ST_prot_anno %>%  mutate(Phenotype = "Serum IL6ST level")
TYK2_prot_anno <- TYK2_prot_anno %>%  mutate(Phenotype = "Serum TYK2 level")
IL6R_prot_anno <- IL6R_prot_anno %>%  mutate(Phenotype = "Serum IL6R level")
IL2RA_prot_anno <- IL2RA_prot_anno %>%  mutate(Phenotype = "Serum IL2RA level")


# Export filest to a CSV file with default parameters of
# TwoSampleMR::read_exposure_data()

#Export
export(IL6ST_prot_anno, file = "data/export_harmonization/IL6ST_prot_anno_TwoSampleMR.csv")
export(IL2RA_prot_anno, file = "data/export_harmonization/IL2RA_prot_anno_TwoSampleMR.csv")
export(IL6R_prot_anno, file = "data/export_harmonization/IL6R_prot_anno_TwoSampleMR.csv")
export(TYK2_prot_anno, file = "data/export_harmonization/TYK2_prot_anno_TwoSampleMR.csv")


#2. eQTL data---------------------------------------------------------------

## 2.1 Import data--------------------------------------------------------------
pacman::p_load(data.table, dtplyr)

names_eqtl <- c("eqtl_IL2RA", "eqtl_IL2RB", "eqtl_IL6R",
                "eqtl_JAK2", "eqtl_JAK3", "eqtl_TYK2", "eqtl_IL6ST")

#Read eqtl sumstat
filenames_eqtl <- paste0("data/export_cis_sumstats/", names_eqtl,"-Gather.rds")
eqtl_sumstats <- map(filenames_eqtl, ~as_tibble(import(file = .x)))

names_eqtl <- c("IL2RA_eqtl", "IL2RB_eqtl", "IL6R_eqtl",
                "JAK2_eqtl",  "JAK3_eqtl", "TYK2_eqtl", "IL6ST_eqtl")

names(eqtl_sumstats) <- names_eqtl

# According to the readme files of eQTLGEN phase I
# (data/REAMDE_cis_eqtl.txt and data/README_cis_eqtl_allele_frequency.txt):
c("Pvalue", "P-value",  "SNP", "SNP rs ID",  "SNPChr", "SNP chromosome",
  "SNPPos", "SNP position",  "AssessedAllele",
  "Assessed allele, the Z-score refers to this allele",  "OtherAllele",
  "Not assessed allele",  "Zscore", "Z-score", "Gene",
  "ENSG name (Ensembl v71) of the eQTL gene",  "GeneSymbol",
  "HGNC name of the gene",  "GeneChr", "Gene chromosome", "GenePos",
  "Centre of gene position",  "NrCohorts",
  "Total number of cohorts where this SNP-gene combination was tested",
  "NrSamples", "Total number of samples where this SNP-gene combination was tested",
  "FDR", "False discovery rate estimated based on permutations",  "BonferroniP",
  "P-value after Bonferroni correction",  "SNP", "SNP rs ID", "hg19_chr",
  "chr number",  "hg19_pos", "SNP position (hg19)",  "AlleleA", "Other allele",
  "AlleleB", "Assessed allele",  "allA_total", "Total allele count of genotype AA",
  "allAB_total", "Total allele count of genotype AB",  "allB_total",
  "Total allele count of genotype BB",  "AlleleB_all",
  "Allele frequency of assessed allele") %>%
  matrix(byrow =T, ncol = 2) %>%
  cbind(c(rep("README_cis_eqtl.txt", 15), rep("README_cis_eqtl_AF.txt", 9))) %>%
  knitr::kable(caption = "sources: README_cis_eqtl_AF.txt = https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_allele_frequency,
  README_cis_eqtl.txt = https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_cis")


# There is some disagreement between the readme and the original file.
# E.g. there are 1,053 subjects near IL2RB which for which 'AssessedAllele'
# is not equal to 'AlleleB'.
colnames(eqtl_sumstats$IL2RB_eqtl)


table_conflicts_eqtl <- c(
  "IL2RA", get_n_conflicts(eqtl_sumstats$IL2RA_eqtl),
  "IL2RB", get_n_conflicts(eqtl_sumstats$IL2RB_eqtl),
  "IL6R", get_n_conflicts(eqtl_sumstats$IL6R_eqtl),
  "JAK2", get_n_conflicts(eqtl_sumstats$JAK2_eqtl),
  "JAK3", get_n_conflicts(eqtl_sumstats$JAK3_eqtl),
  "TYK2", get_n_conflicts(eqtl_sumstats$TYK2_eqtl),
  "IL6ST", get_n_conflicts(eqtl_sumstats$TYK2_eqtl)) %>%
  matrix(ncol = 2, byrow = TRUE)

colnames(table_conflicts_eqtl) <- c("eQTL", "n conflicts")

knitr::kable(table_conflicts_eqtl)



# In a remarkably stable pattern between files, AlleleB == Assessed allele when
# whenAF < 0.50 and AlleleB != Assessed allele when AF > 0.50.
map(eqtl_sumstats, summarize_agreement)

# From this, it is taken that AF is given per effect alleles mentioned
# in README_cis_eqtl_AF.txt. Therefore, EAFs are harmonized to align with
# the variable 'AssessedAllele'
eqtl_sumstats <- map(eqtl_sumstats, ~.x %>% rename(AF_B = AlleleB_all))
eqtl_sumstats <- map(eqtl_sumstats, ~.x %>%
                       mutate(EAF =
                                ifelse(AlleleB == AssessedAllele,
                                       AF_B, 1 - AF_B)))
eqtl_sumstats <- map(eqtl_sumstats, ~.x %>%
  mutate(MAF = ifelse(EAF < 0.5,  EAF, 1 - EAF)))

#Calculate betas and SEs
eqtl_sumstats <- map(eqtl_sumstats, ~.x %>%
                       mutate(beta = Zscore / sqrt(2 * AF_B * (1-AF_B) * (NrSamples + Zscore^2)),
                              se = 1 / sqrt(2 * AF_B * (1-AF_B) * (NrSamples + Zscore^2))))
#This way of calculating betas and SEs is mentioned in summary data –based MR-section #(https://www.eqtlgen.org/cis-eqtls.html) README file

####"Calculations
#### #------------#
#### The minor allele frequencies were calculated using reported allele counts from all cohorts except
#### Framingham Heart Study, because the relatedness in this cohort may influence allele frequencies.
#### Beta and standard error values were calculated as follows:
#### beta = Z-score / sqrt(2 * alleleFreq * (1-alleleFreq) * (N + Z-score^2)
#### se = 1 / sqrt(2 * alleleFreq * (1-alleleFreq) * (N + Z-score^2)
####
#### Where N = sample size (max. 31864)"

#Rename files
eqtl_traitnames <- paste0("Whole blood ",
                          c("IL2RA", "IL2RB", "IL6R", "JAK2", "JAK3", "TYK2", "IL6ST"),
                          " gene expression")

eqtl_sumstats <- map2(eqtl_sumstats, eqtl_traitnames, ~.x %>% mutate(Phenotype = .y))
eqtl_sumstats <- map(eqtl_sumstats, rename_eqtl_data)

#Export to csv-files
eqtl_exportnames <- paste0("data/export_harmonization/", names_eqtl, "_TwoSampleMR.csv")
map2(eqtl_sumstats, eqtl_exportnames, ~export(.x, file = .y))


#3. Outcome data---------------------------------------------------------------

#3.1 Import outcome data-------------------------------------------------------
TYK2_T1D    <- import("data/export_cis_sumstats/TYK2_T1D.rds")  %>% as_tibble()
JAK2_T1D    <- import("data/export_cis_sumstats/JAK2_T1D.rds")  %>% as_tibble()
JAK3_T1D    <- import("data/export_cis_sumstats/JAK3_T1D.rds")  %>% as_tibble()
IL6ST_T1D   <- import("data/export_cis_sumstats/IL6ST_T1D.rds") %>% as_tibble()
IL6R_T1D    <- import("data/export_cis_sumstats/IL6R_T1D.rds")  %>% as_tibble()
IL2RB_T1D   <- import("data/export_cis_sumstats/IL2RB_T1D.rds") %>% as_tibble()
IL2RA_T1D   <- import("data/export_cis_sumstats/IL2RA_T1D.rds") %>% as_tibble()

# There were no missing rsid values
TYK2_T1D   %>% pull(hm_rsid) %>% anyNA()
JAK2_T1D   %>% pull(hm_rsid) %>% anyNA()
JAK3_T1D   %>% pull(hm_rsid) %>% anyNA()
IL6ST_T1D  %>% pull(hm_rsid) %>% anyNA()
IL6R_T1D   %>% pull(hm_rsid) %>% anyNA()
IL2RB_T1D  %>% pull(hm_rsid) %>% anyNA()
IL2RA_T1D  %>% pull(hm_rsid) %>% anyNA()


## 3.2 Remove duplicated SNPs--------------------------------------------------

# All data frames had duplicated rsids (suppressed). They are
# removed.
TYK2_T1D   %>% pull(hm_rsid) %>% anyDuplicated()
JAK2_T1D   %>% pull(hm_rsid) %>% anyDuplicated()
JAK3_T1D   %>% pull(hm_rsid) %>% anyDuplicated()
IL6ST_T1D  %>% pull(hm_rsid) %>% anyDuplicated()
IL6R_T1D   %>% pull(hm_rsid) %>% anyDuplicated()
IL2RB_T1D  %>% pull(hm_rsid) %>% anyDuplicated()
IL2RA_T1D  %>% pull(hm_rsid) %>% anyDuplicated()

duplicated_snps_TYK2   <- TYK2_T1D   %>%
  filter(duplicated(hm_rsid)) %>% pull(hm_rsid)
duplicated_snps_JAK2   <- JAK2_T1D   %>%
  filter(duplicated(hm_rsid)) %>% pull(hm_rsid)
duplicated_snps_JAK3   <- JAK3_T1D   %>%
  filter(duplicated(hm_rsid)) %>% pull(hm_rsid)
duplicated_snps_IL6ST  <- IL6ST_T1D  %>%
  filter(duplicated(hm_rsid)) %>% pull(hm_rsid)
duplicated_snps_IL6R   <- IL6R_T1D   %>%
  filter(duplicated(hm_rsid)) %>% pull(hm_rsid)
duplicated_snps_IL2RB  <- IL2RB_T1D  %>%
  filter(duplicated(hm_rsid)) %>% pull(hm_rsid)
duplicated_snps_IL2RA  <- IL2RA_T1D  %>%
  filter(duplicated(hm_rsid)) %>% pull(hm_rsid)

TYK2_T1D   <- TYK2_T1D    %>% filter(!hm_rsid %in% duplicated_snps_TYK2)
JAK2_T1D   <- JAK2_T1D    %>% filter(!hm_rsid %in% duplicated_snps_JAK2)
JAK3_T1D   <- JAK3_T1D    %>% filter(!hm_rsid %in% duplicated_snps_JAK3)
IL6ST_T1D  <- IL6ST_T1D %>% filter(!hm_rsid %in% duplicated_snps_IL6ST)
IL6R_T1D   <- IL6R_T1D    %>% filter(!hm_rsid %in% duplicated_snps_IL6R)
IL2RB_T1D  <- IL2RB_T1D   %>% filter(!hm_rsid %in% duplicated_snps_IL2RB)
IL2RA_T1D  <- IL2RA_T1D   %>% filter(!hm_rsid %in% duplicated_snps_IL2RA)


## 3.3 Write data in TwoSampleMR format-----------------------------------------

#Change to numeric


IL2RA_T1D   <- IL2RA_T1D  %>% rename_for_TwoSampleMR_o()
IL2RB_T1D   <- IL2RB_T1D  %>% rename_for_TwoSampleMR_o()
IL6ST_T1D   <- IL6ST_T1D  %>% rename_for_TwoSampleMR_o()
IL6R_T1D    <- IL6R_T1D   %>% rename_for_TwoSampleMR_o()
JAK2_T1D    <- JAK2_T1D   %>% rename_for_TwoSampleMR_o()
JAK3_T1D    <- JAK3_T1D   %>% rename_for_TwoSampleMR_o()
TYK2_T1D    <- TYK2_T1D   %>% rename_for_TwoSampleMR_o()

# Write data in TwoSampleMR data format
IL2RA_T1D  %>%
  export("data/export_harmonization/IL2RA_T1D_TwoSampleMR.csv")
IL2RB_T1D  %>%
  export("data/export_harmonization/IL2RB_T1D_TwoSampleMR.csv")
IL6ST_T1D  %>%
  export("data/export_harmonization/IL6ST_T1D_TwoSampleMR.csv")
IL6R_T1D  %>%
  export("data/export_harmonization/IL6R_T1D_TwoSampleMR.csv")
JAK2_T1D   %>%
  export("data/export_harmonization/JAK2_T1D_TwoSampleMR.csv")
JAK3_T1D   %>%
  export("data/export_harmonization/JAK3_T1D_TwoSampleMR.csv")
TYK2_T1D   %>%
  export("data/export_harmonization/TYK2_T1D_TwoSampleMR.csv")

#Sanity check: no warnings
# read_outcome_data("data/export_harmonization/IL2RA_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/IL2RB_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/IL6ST_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/IL6R_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/JAK2_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/JAK3_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/TYK2_T1D_TwoSampleMR.csv", sep = ",")


#4 Tissue atlas eQTL--------------------------------------------------------------

IL2RA_spleen <-
  read_tsv(here::here("data/import/singlecellqtl-VK-May3/IL2RA_spleen.tsv")) %>%
  mutate(tissue = "Spleen", gene = "IL2RA")
IL6R_spleen <-
  read_tsv(here::here("data/import/singlecellqtl-VK-May3/IL6R_spleen.tsv")) %>%
  mutate(tissue = "Spleen", gene = "IL6R")
IL6R_pancreas <-
  read_tsv(here::here("data/import/singlecellqtl-VK-May3/IL6R_pancreas.tsv")) %>%
  mutate(tissue = "Pancreas", gene = "IL6R")
TYK2_spleen <-
  read_tsv(here::here("data/import/singlecellqtl-VK-May3/TYK2_spleen.tsv")) %>%
  mutate(tissue = "Spleen", gene = "TYK2")
TYK2_pancreas <-
  read_tsv(here::here("data/import/singlecellqtl-VK-May3/TYK2_pancreas.tsv")) %>%
  mutate(tissue = "Pancreas", gene = "TYK2")

IL6ST_spleen <-
  read_tsv(here::here("data/import/singlecellqtl-VK-May3/IL6ST_spleen_noheader.tsv"),
           col_names = FALSE) %>%
  mutate(tissue = "Spleen", gene = "IL6ST") %>%
  set_names(colnames(TYK2_pancreas))
IL6ST_pancreas <-
  read_tsv(here::here("data/import/singlecellqtl-VK-May3/IL6ST_pancreas_noheader.tsv"),
           col_names = FALSE) %>%
  mutate(tissue = "Pancreas", gene = "IL6ST")  %>%
  set_names(colnames(TYK2_pancreas))

atlas <- IL2RA_spleen %>%
  bind_rows(IL6R_spleen) %>%
  bind_rows(IL6R_pancreas) %>%
  bind_rows(TYK2_spleen) %>%
  bind_rows(TYK2_pancreas) %>%
  bind_rows(IL6ST_spleen) %>%
  bind_rows(IL6ST_pancreas)

rm(IL2RA_spleen, IL6R_spleen, IL6R_pancreas, TYK2_spleen, TYK2_pancreas,
   IL6ST_pancreas, IL6ST_spleen)

atlas <- atlas %>%
  mutate(samplesize = an/2) %>%
  rename(pos = position, effect_allele = ref, other_allele = alt,
         SNP = rsid, maf = maf, pval = pvalue) %>%
  mutate(Phenotype = paste0(tissue, " ", gene, " mRNA expression")) %>%
  filter(!is.na(SNP))#remove 8 unname d SNPs

#Number of samples
atlas %>%
  distinct(gene, tissue, samplesize) %>%
  group_by(gene, tissue) %>%
  summarise(min_samplesize = min(samplesize),
            max_samplesize = max(samplesize))

#Export data
atlas  %>%
  group_by(tissue, gene) %>%
  group_split() %>%
  iwalk(~write_csv(.x, paste0(here("data/export_celltypes_harmonized/"), .x$gene[1], "_", .x$tissue[1], "_TwoSampleMR.csv")))


# exposure_dat <- read_exposure_data("data/export_harmonization/IL2RA_spleen_TwoSampleMR.csv", sep = ",")

#5 Single-cell eQTL-------------------------------------------------------------


#1. Read orignal files and combine them

library(data.table)
# eqtl_table <- read_tsv(here::here("data/import/eqtl_table.tsv"))
# eqtl_IL6ST <- eqtl_table %>%  filter(gene == "IL6ST")

#The previous command would work, but is less memory efficient than
# awk -F'\t' '$5 == "IL6ST"' eqtl_table.tsv > eqtl_table_IL6ST.tsv

# awk -F'\t' 'NR == 1 || $5 == "IL6ST" || $5 == "IL2RA" || $5 == "TYK2" || $5 == "IL6R"' eqtl_table.tsv > eqtl_table_targets.tsv


IL6ST_sc <- read_tsv(here::here("data/import/eqtl_table_IL6ST.tsv"))

#Read
IL2RA_sc <- read_tsv(
  here::here(
    "data/import/singlecellqtl-JL-Apr29/input_scQTL/IL2RA_single_cellQTL.tsv"))
IL6R_sc <- read_tsv(
  here::here(
    "data/import/singlecellqtl-JL-Apr29/input_scQTL/IL6R_single_cellQTL.tsv"))
TYK2_sc <- read_tsv(
  here::here(
    "data/import/singlecellqtl-JL-Apr29/input_scQTL/TYK2_single_cellQTL.tsv"))

#Combine
sc <- IL2RA_sc %>%
  bind_rows(IL6R_sc) %>%
  bind_rows(TYK2_sc) %>%
  bind_rows(IL6ST_sc)

#Add N of each cell type / participant based on Jake Lin's email on Apr 29, 2024:
sc_N <- sc %>%
  distinct(CELL_ID, CELL_TYPE) %>%
  mutate(samplesize = case_when(
    CELL_ID == "cd8s100b" ~ 981,
    CELL_ID == "monoc" ~ 969,
    CELL_ID == "nk" ~ 928,
    CELL_ID == "nkr" ~ 967,
    CELL_ID == "mononc" ~ 932,
    CELL_ID == "plasma" ~ 643,
    T ~ 982))
sc <- left_join(sc, sc_N, by = c("CELL_ID", "CELL_TYPE"))

#Remove originals
rm(IL2RA_sc, IL6R_sc, TYK2_sc, IL6ST_sc)

#Get the number of variants for each combination
sc %>%
  filter(!is.na(P_VALUE)) %>%
  group_by(GENE, CELL_ID, CELL_TYPE) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = GENE, values_from = n)

#2. Rename data.frames for read_exposure_data() and save as csv
sc %>%
  rename(pos = POS, effect_allele = A2, other_allele = A1,
         SNP = RSID, eaf = A2_FREQ_ONEK1K, pval = P_VALUE) %>%
  mutate(Phenotype = paste0(CELL_TYPE, " ", GENE, " mRNA expression")) %>%
  group_by(CELL_ID, GENE) %>%
  group_split() %>%
  iwalk(~write_csv(.x, paste0(here("data/export_celltypes_harmonized/"), .x$GENE[1], "_", .x$CELL_ID[1], "_TwoSampleMR.csv")))

sc_names <-  sc %>%  distinct(CELL_TYPE, CELL_ID) %>% rename(celltype = CELL_ID, celltype_long = CELL_TYPE)
write_csv(sc_names, here::here("data/export_harmonization/single_cell_types.csv"))

