#Load package



# This script merges the annotation and protein data.
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

pqtl_names <- c("IL2RA", "IL2RB", "IL2RG", "IL6R",
                "IL6ST", "IL12B", "JAK2",  "TYK2")

#Load GWAS for proteins
pqtl_prot_names <- c("IL2RA_prot", "IL2RB_prot", "IL2RG_prot", "IL6R_prot",
                "IL6ST_prot", "IL12B_prot", "JAK2_prot",  "TYK2_prot")
pqtl_prot_filenames <- paste0("data/cis_sumstats/", pqtl_prot_names, ".rds")
pqtl_prot <- map(pqtl_prot_filenames, ~import(.x))

#Load further annotation for GWAS of pqtl
pqtl_anno_names <- c("IL2RA_anno", "IL2RB_anno", "IL2RG_anno", "IL6R_anno",
                "IL6ST_anno", "IL12B_anno", "JAK2_anno",  "TYK2_anno")
pqtl_anno_filenames <- paste0("data/cis_sumstats/", pqtl_anno_names, ".rds")
pqtl_anno <- map(pqtl_anno_filenames, ~import(.x))


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

left_join_pqtl <- function(x, y) left_join(x, y, by = c("Name", "Chrom", "Pos"))
pqtl <- map2(pqtl_prot, pqtl_anno, left_join_pqtl)
names(pqtl) <- pqtl_names

## 1.3 Sort conflicts-----------------------------------------------------------

#Further in the README file of Ferkingstad et al.
#
#### Also, the summary statistics files sometimes incorrectly have effectAllele=otherAllele for multiallelic variants.
#### In these cases the effectAllele is correct, but the otherAllele should be '!', meaning that the effectAllele is tested against the other (two or more) alleles (using the '!' sign as shorthand for 'not effectAllele').
#### This has been corrected in the file Extra annotation file (assocvariants.annotated.txt.gz).
#
#Therefore, multiallelic variants were removed.

filter_multiallelic <- function(x) filter(x, otherAllele.y != "!")
pqtl <- map(pqtl, filter_multiallelic)


# After removal of multiallelic variants there were no conflicts
# between effectAllele.x and effectAllele.y or between
# otherAllele.x and otherAllele.y

#Following function is simply to avoid repetition
check_conflicts_ea_oa <- function(data, effectAllele.x, effectAllele.y) {
  filter(data, effectAllele.x != effectAllele.y |
           is.na(effectAllele.x) & !is.na(effectAllele.y) |
           !is.na(effectAllele.x) &  is.na(effectAllele.y) |
           otherAllele.x != otherAllele.y |
            is.na(otherAllele.x) & !is.na(otherAllele.y) |
           !is.na(otherAllele.x) &  is.na(otherAllele.y))
}

#List those with conflicts in effectAllele or otherAllele
map(pqtl, check_conflicts_ea_oa) #none

# After this, the remaining discrepancies between the variables
# obtained from annotation and protein data files were because those
# variants without an rsid were listed as "." in the annotation file
# and as "NA" in the main GWAS sumstat file

#List subjects which had conflicting rsid information and show
#which entries there are for these
check_conflicts_rsid <- function(data, effectAllele.x, effectAllele.y) {
  filter(data, rsids.x != rsids.y |
           is.na(rsids.x) & !is.na(rsids.y) |
           !is.na(rsids.x) &  is.na(rsids.y)) %>%
    mutate(rsids.x = factor(rsids.x),
         rsids.y = factor(rsids.y)) %>%
  select(rsids.x, rsids.y) %>%
  summary()
}

#All the conflicts in rsid.x and rsid.y are because rsid.x is NA and rsid.y is "."
map(pqtl, check_conflicts_rsid)

# Thus, only those values variables with .x in their name are preserved.
# Here, we also remove the remaining multiallelic SNPs
# (where effect allele equals the other allele.)


#Discard all variables that end with an ".y"
pqtl <- map(pqtl, ~.x %>% select(!ends_with(".y")))

#Remove ".x" from variable names
pqtl <- map(pqtl, ~.x %>% rename_with(~gsub(".x", "", .x, fixed = TRUE)))

#Sanity checks
map(pqtl, colnames)


# All data files seem to have entries without rsids. They are now removed
map(pqtl, ~anyNA(.x$rsids))
pqtl <- map(pqtl, ~filter(.x, !is.na(rsids)))
# map(pqtl, ~anyNA(.x$rsids)) #none

## 1.4 Write data in TwoSampleMR format-----------------------------------------

# Data are written in TwoSampleMR format so that they can be easily opened
# in the main analyses scripts.

# Rename files
rename_for_TwoSampleMR_e  <- function(data, rsids, Beta, SE, effectAlleleFreq,
                                    effectAllele, otherAllele,
                                    Pval, Pos, N) {
  data %>%
    rename(SNP = rsids, beta = Beta, se = SE, eaf = effectAlleleFreq,
         effect_allele = effectAllele, other_allele = otherAllele,
         pval = Pval, pos = Pos) %>%
    mutate(samplesize = N)
  }

#Rename
pqtl <- map(pqtl, rename_for_TwoSampleMR_e)

#Add phenotype information
pqtl_phenotype_names <- paste0("Serum ", pqtl_names, " level")
pqtl <- map2(pqtl, pqtl_phenotype_names, ~.x %>% mutate(Phenotype = .y))


# Export filest to a CSV file with default parameters of
# TwoSampleMR::read_exposure_data()

#Export
pqtl_filenames <-  paste0("data/export_harmonization/", pqtl_names, "_prot_anno_TwoSampleMR.csv")
map2(pqtl, pqtl_filenames, ~export(x=.x, file = .y))

#Sanity check: no warnings
# map(pqtl_filenames, ~read_exposure_data("data/export_harmonization/IL12B_prot_anno_TwoSampleMR.csv", sep = ","))


#2. IL6 CRP data---------------------------------------------------------------

## 2.1 Write IL6R_crp in TwoSampleMR format-------------------------------------

IL6R_crp   <- import("data/cis_sumstats/IL6R_crp.rds")

# colnames(IL6R_crp)
IL6R_crp <- IL6R_crp %>%
  select(-other_allele, -effect_allele, -beta) %>%
  rename(SNP = hm_rsid, beta = hm_beta, se = standard_error,
           eaf = hm_effect_allele_frequency,
           effect_allele = hm_effect_allele,
           other_allele = hm_other_allele,
           pval = p_value, pos = hm_pos) %>%
  mutate(pval = as.numeric(pval),
         eaf = as.numeric(eaf),
         samplesize = 575531,
         Phenotype = "Serum C-reactive protein levels")

IL6R_crp  %>%
  export("data/export_harmonization/IL6R_crp_TwoSampleMR.csv")

#eQTL data

## 2.2 Import data--------------------------------------------------------------
pacman::p_load(data.table, dtplyr)

names_eqtl <- c("eqtl_IFNAR2", "eqtl_IL2RA", "eqtl_IL2RG", "eqtl_IL2RB",
                "eqtl_IL6R", "eqtl_JAK1", "eqtl_JAK2", "eqtl_JAK3", "eqtl_TYK2")

#Read eqtl sumstat
filenames_eqtl <- paste0("data/cis_eqtl_sumstats/", names_eqtl,"-Gather.rds")
eqtl_sumstats <- map(filenames_eqtl, ~as_tibble(import(file = .x)))

names_eqtl <- c("IFNAR2_eqtl",  "IL2RA_eqtl", "IL2RG_eqtl", "IL2RB_eqtl",
                "IL6R_eqtl",  "JAK1_eqtl", "JAK2_eqtl",  "JAK3_eqtl",
                "TYK2_eqtl")

names(eqtl_sumstats) <- names_eqtl
#IFNAR2 had a lot of duplicated SNPs, and on closer inspection, there are two genes with different Ensmbl-symbols under IFNAR2 (ENSG00000159110 and ENSG00000249624). The former is IFNAR2 and the other codes a "novel gene"
eqtl_sumstats$IFNAR2_eqtl$Gene %>% factor() %>% levels()

#None of the other genes had two Ensemble codes (not shown)

#Therefore, all the entries with "Gene == ENSG00000249624" were removed
eqtl_sumstats$IFNAR2_eqtl <- eqtl_sumstats$IFNAR2_eqtl %>% filter(Gene == "ENSG00000159110")

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

get_n_conflicts  <- function(data, AlleleB=AlleleB, AssessedAllele=AssessedAllele) {
  data %>%
  filter(AlleleB != AssessedAllele) %>%
  pull(SNP) %>% length()
}

table_conflicts_eqtl <- c("IFNAR2", get_n_conflicts(eqtl_sumstats$IFNAR2_eqtl),
  "IL2RA", get_n_conflicts(eqtl_sumstats$IL2RA_eqtl),
  "IL2RG", get_n_conflicts(eqtl_sumstats$IL2RG_eqtl),
  "IL2RB", get_n_conflicts(eqtl_sumstats$IL2RB_eqtl),
  "IL6R", get_n_conflicts(eqtl_sumstats$IL6R_eqtl),
  "JAK1", get_n_conflicts(eqtl_sumstats$JAK1_eqtl),
  "JAK2", get_n_conflicts(eqtl_sumstats$JAK2_eqtl),
  "JAK3", get_n_conflicts(eqtl_sumstats$JAK3_eqtl),
  "TYK2", get_n_conflicts(eqtl_sumstats$TYK2_eqtl)) %>%
  matrix(ncol = 2, byrow = TRUE)

colnames(table_conflicts_eqtl) <- c("eQTL", "n conflicts")

knitr::kable(table_conflicts_eqtl)


summarize_agreement <- function(data, AlleleB=AlleleB,
                           AssessedAllele=AssessedAllele,
                           AlleleB_all=AlleleB_all) {
  data %>%
  mutate(AlleleB_vs_AssessedAllele =
           ifelse(AlleleB == AssessedAllele, "AlleleB == AssessedAllele",
                  "AlleleB != AssessedAllele")) %>%
  group_by(AlleleB_vs_AssessedAllele) %>%
  summarise(p_AF_m40  = round(sum(AlleleB_all < 0.40) / length(AlleleB_all), 4),
            p_AF_4050 = round(sum(AlleleB_all >= 0.40 & AlleleB_all < 0.50)  /
                                length(AlleleB_all), 4),
            p_AF_5060 = round(sum(AlleleB_all > 0.50 & AlleleB_all <= 0.60)  /
                                length(AlleleB_all), 4),
            p_AF_p60  = round(sum(AlleleB_all > 0.60) / length(AlleleB_all), 4),
            p_AF_50  = round(sum(AlleleB_all == 0.50) / length(AlleleB_all), 4),
            n_AF_NA = sum(is.na(AlleleB_all)))
}

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
eqtl_traitnames <- paste0("Whole blood ", c("IFNAR2",  "IL2RA",
                     "IL2RG", "IL2RB",
                     "IL6R", "JAK1",
                     "JAK2", "JAK3",
                     "TYK2"), " gene expression")

rename_eqtl_data <- function(data, Phenotype, EAF, AssessedAllele, OtherAllele,
                              Pvalue, SNPPos, NrSamples) {
   data %>%
  rename(eaf = EAF, effect_allele = AssessedAllele,
         other_allele = OtherAllele, pval = Pvalue,
         pos = SNPPos, chr = SNPChr, samplesize = NrSamples) %>%
  mutate(pval = as.numeric(pval),
         eaf = as.numeric(eaf),
         Phenotype = Phenotype)
 }
eqtl_sumstats <- map2(eqtl_sumstats, eqtl_traitnames, ~.x %>% mutate(Phenotype = .y))
eqtl_sumstats <- map(eqtl_sumstats, rename_eqtl_data)

#Export to csv-files
eqtl_exportnames <- paste0("data/export_harmonization/", names_eqtl, "_TwoSampleMR.csv")
map2(eqtl_sumstats, eqtl_exportnames, ~export(.x, file = .y))


#3. Outcome data---------------------------------------------------------------

#3.1 Import outcome data-------------------------------------------------------
IL12B_T1D   <- import("data/cis_sumstats/IL12B_T1D.rds") %>% as_tibble()
TYK2_T1D    <- import("data/cis_sumstats/TYK2_T1D.rds")  %>% as_tibble()
JAK1_T1D    <- import("data/cis_sumstats/JAK1_T1D.rds")  %>% as_tibble()
JAK2_T1D    <- import("data/cis_sumstats/JAK2_T1D.rds")  %>% as_tibble()
JAK3_T1D    <- import("data/cis_sumstats/JAK3_T1D.rds")  %>% as_tibble()
IL6ST_T1D   <- import("data/cis_sumstats/IL6ST_T1D.rds") %>% as_tibble()
IL6R_T1D    <- import("data/cis_sumstats/IL6R_T1D.rds")  %>% as_tibble()
IL2RB_T1D   <- import("data/cis_sumstats/IL2RB_T1D.rds") %>% as_tibble()
IL2RG_T1D   <- import("data/cis_sumstats/IL2RG_T1D.rds") %>% as_tibble()
IL2RA_T1D   <- import("data/cis_sumstats/IL2RA_T1D.rds") %>% as_tibble()
IFNAR2_T1D  <- import("data/cis_sumstats/IFNAR2_T1D.rds") %>% as_tibble()

# There were no missing rsid values (suppressed)

#every
IL12B_T1D  %>% pull(hm_rsid) %>% anyNA()
TYK2_T1D   %>% pull(hm_rsid) %>% anyNA()
JAK1_T1D   %>% pull(hm_rsid) %>% anyNA()
JAK2_T1D   %>% pull(hm_rsid) %>% anyNA()
JAK3_T1D   %>% pull(hm_rsid) %>% anyNA()
IL6ST_T1D  %>% pull(hm_rsid) %>% anyNA()
IL6R_T1D   %>% pull(hm_rsid) %>% anyNA()
IL2RB_T1D  %>% pull(hm_rsid) %>% anyNA()
IL2RG_T1D  %>% pull(hm_rsid) %>% anyNA()
IL2RA_T1D  %>% pull(hm_rsid) %>% anyNA()
IFNAR2_T1D %>% pull(hm_rsid) %>% anyNA()


## 3.2 Remove duplicated SNPs--------------------------------------------------

# All data frames had duplicated rsids (suppressed). They are
# removed.
IL12B_T1D  %>% pull(hm_rsid) %>% anyDuplicated()
TYK2_T1D   %>% pull(hm_rsid) %>% anyDuplicated()
JAK1_T1D   %>% pull(hm_rsid) %>% anyDuplicated()
JAK2_T1D   %>% pull(hm_rsid) %>% anyDuplicated()
JAK3_T1D   %>% pull(hm_rsid) %>% anyDuplicated()
IL6ST_T1D  %>% pull(hm_rsid) %>% anyDuplicated()
IL6R_T1D   %>% pull(hm_rsid) %>% anyDuplicated()
IL2RB_T1D  %>% pull(hm_rsid) %>% anyDuplicated()
IL2RG_T1D  %>% pull(hm_rsid) %>% anyDuplicated()
IL2RA_T1D  %>% pull(hm_rsid) %>% anyDuplicated()
IFNAR2_T1D %>% pull(hm_rsid) %>% anyDuplicated()

duplicated_snps_IL12B  <- IL12B_T1D  %>%
  filter(duplicated(hm_rsid)) %>% pull(hm_rsid)
duplicated_snps_TYK2   <- TYK2_T1D   %>%
  filter(duplicated(hm_rsid)) %>% pull(hm_rsid)
duplicated_snps_JAK1   <- JAK1_T1D   %>%
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
duplicated_snps_IL2RG  <- IL2RG_T1D  %>%
  filter(duplicated(hm_rsid)) %>% pull(hm_rsid)
duplicated_snps_IL2RA  <- IL2RA_T1D  %>%
  filter(duplicated(hm_rsid)) %>% pull(hm_rsid)
duplicated_snps_IFNAR2 <- IFNAR2_T1D  %>%
  filter(duplicated(hm_rsid)) %>% pull(hm_rsid)

IL12B_T1D  <- IL12B_T1D   %>% filter(!hm_rsid %in% duplicated_snps_IL12B)
TYK2_T1D   <- TYK2_T1D    %>% filter(!hm_rsid %in% duplicated_snps_TYK2)
JAK1_T1D   <- JAK1_T1D    %>% filter(!hm_rsid %in% duplicated_snps_JAK1)
JAK2_T1D   <- JAK2_T1D    %>% filter(!hm_rsid %in% duplicated_snps_JAK2)
JAK3_T1D   <- JAK3_T1D    %>% filter(!hm_rsid %in% duplicated_snps_JAK3)
IL6ST_T1D  <- IL6ST_T1D %>% filter(!hm_rsid %in% duplicated_snps_IL6ST)
IL6R_T1D   <- IL6R_T1D    %>% filter(!hm_rsid %in% duplicated_snps_IL6R)
IL2RB_T1D  <- IL2RB_T1D   %>% filter(!hm_rsid %in% duplicated_snps_IL2RB)
IL2RG_T1D  <- IL2RG_T1D   %>% filter(!hm_rsid %in% duplicated_snps_IL2RG)
IL2RA_T1D  <- IL2RA_T1D   %>% filter(!hm_rsid %in% duplicated_snps_IL2RA)
IFNAR2_T1D <- IFNAR2_T1D  %>% filter(!hm_rsid %in% duplicated_snps_IFNAR2)


## 3.3 Write data in TwoSampleMR format-----------------------------------------

#Change to numeric
rename_for_TwoSampleMR_o  <- function(data, beta, other_allele, effect_allele,
                                      hm_rsid, hm_beta, standard_error,
                                      hm_effect_allele_frequency,
                                      hm_effectAllele, hm_otherAllele,
                                      p_value, hm_pos) {
  data %>%
  select(-beta, -other_allele, -effect_allele) %>%
  rename(SNP = hm_rsid, beta = hm_beta, se = standard_error, eaf = hm_effect_allele_frequency,
         effect_allele = hm_effect_allele, other_allele = hm_other_allele,
         pval = p_value, pos = hm_pos) %>%
  mutate(Phenotype = "Risk of type 1 diabetes",
         ncase = 18942, ncontrol = 501638, samplesize = 18942+501638)
  }

IL12B_T1D   <- IL12B_T1D  %>% rename_for_TwoSampleMR_o()
IL2RA_T1D   <- IL2RA_T1D  %>% rename_for_TwoSampleMR_o()
IL2RB_T1D   <- IL2RB_T1D  %>% rename_for_TwoSampleMR_o()
IL2RG_T1D   <- IL2RG_T1D  %>% rename_for_TwoSampleMR_o()
IL6ST_T1D   <- IL6ST_T1D  %>% rename_for_TwoSampleMR_o()
IL6R_T1D    <- IL6R_T1D   %>% rename_for_TwoSampleMR_o()
JAK1_T1D    <- JAK1_T1D   %>% rename_for_TwoSampleMR_o()
JAK2_T1D    <- JAK2_T1D   %>% rename_for_TwoSampleMR_o()
JAK3_T1D    <- JAK3_T1D   %>% rename_for_TwoSampleMR_o()
TYK2_T1D    <- TYK2_T1D   %>% rename_for_TwoSampleMR_o()
IFNAR2_T1D  <- IFNAR2_T1D %>% rename_for_TwoSampleMR_o()



# Write data in TwoSampleMR data format
IL12B_T1D  %>%
  export("data/export_harmonization/IL12B_T1D_TwoSampleMR.csv")
IL2RA_T1D  %>%
  export("data/export_harmonization/IL2RA_T1D_TwoSampleMR.csv")
IL2RB_T1D  %>%
  export("data/export_harmonization/IL2RB_T1D_TwoSampleMR.csv")
IL2RG_T1D  %>%
  export("data/export_harmonization/IL2RG_T1D_TwoSampleMR.csv")
IL6ST_T1D  %>%
  export("data/export_harmonization/IL6ST_T1D_TwoSampleMR.csv")
IL6R_T1D  %>%
  export("data/export_harmonization/IL6R_T1D_TwoSampleMR.csv")
JAK1_T1D   %>%
  export("data/export_harmonization/JAK1_T1D_TwoSampleMR.csv")
JAK2_T1D   %>%
  export("data/export_harmonization/JAK2_T1D_TwoSampleMR.csv")
JAK3_T1D   %>%
  export("data/export_harmonization/JAK3_T1D_TwoSampleMR.csv")
TYK2_T1D   %>%
  export("data/export_harmonization/TYK2_T1D_TwoSampleMR.csv")
IFNAR2_T1D   %>%
  export("data/export_harmonization/IFNAR2_T1D_TwoSampleMR.csv")

#Sanity check: no warnings
# read_outcome_data("data/export_harmonization/IL12B_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/IL2RA_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/IL2RB_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/IL2RG_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/IL6ST_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/IL6R_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/JAK1_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/JAK2_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/JAK3_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/TYK2_T1D_TwoSampleMR.csv", sep = ",")
# read_outcome_data("data/export_harmonization/IFNAR2_T1D_TwoSampleMR.csv", sep = ",")


