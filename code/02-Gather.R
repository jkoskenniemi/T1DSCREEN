
mbp <- 1000000

# This script was run in a HPC separately to get the source data.
# setwd("/scratch/project_2006579/") #set working directory in our local HPC

# Changes Apr 25, 2023
# 1. An attempt is made to make the entire analysis workflow traceable.
# 2. Export path is changed from /data/sumstats/ to /data/cis_sumstats/
# 3. Uncommented and modified IL2RG (protein, annotation & T1D sumstats)

#1 Read circulating protein level of IL6ST-----------------------------------------

# Summary statistics of proteins can be downloaded at
# https://www.decode.com/summarydata/ under 2021 and
#
# Ferkingstad, E. et al. Large-scale integration of the plasma proteome with genetics and disease. Summary data (list).
# Extra annotation (430 MB). Excluded variants (52 MB)  Read me file (2 KB).
#
# Direct links to downloads (you have to fill a form first)
# Summary data (list): https://download.decode.is/form/folder/proteomics
#    -> 2620_4_IL6ST_gp130__soluble.txt.gz (decompress to get to 2620_4_IL6ST_gp130__soluble.txt)
#    -> save to data/import/
# Extra annotation: https://download.decode.is/form/2021/assocvariants.annotated.txt.gz
#    -> assocvariants.annotated.txt.gz (decompress to get to assocvariants.annotated.txt)
#    -> save to data/import/
# Readme file: https://download.decode.is/form/2021/proteomics_readme.txt
#    -> proteomics_readme.txt
#    -> save to data/import/

# Read annotations of Ferkingstad et al. 2021 data.

#This file originates from
anno  <- fread("data/import/assocvariants.annotated.txt") #missing

#IL2RA

#Read data
cat("Start parsing IL6ST")

#IL6ST
IL6ST <- fread("data/import/2620_4_IL6ST_gp130__soluble.txt")
IL6ST <- IL6ST[Chrom == "chr5"]
IL6ST <- IL6ST[Pos >=  55230923 - mbp & Pos < 55290821 + mbp]
IL6ST_anno <- anno[Chrom == "chr5"]
IL6ST_anno <- IL6ST_anno[Pos >=  55230923 - mbp & Pos <  55290821 + mbp]
saveRDS(IL6ST, file = "data/cis_sumstats/IL6ST_prot.rds")
saveRDS(IL6ST_anno, file = "data/cis_sumstats/IL6ST_anno.rds")

rm(list=ls())
gc()



  #2 Read data of type 1 diabetes------------------------------------------------

  #Data from GWAS catalog
  #README: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014023/README.txt
  #http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014023/harmonised/34012112-GCST90014023-EFO_0001359.h.tsv.gz
  #   -> decompress and save to data/import/34012112-GCST90014023-EFO_0001359.h.tsv

  T1D <- fread("data/import/34012112-GCST90014023-EFO_0001359.h.tsv")

  #3 Filter areas of interest in GWAS of T1D for each protein--------------------
  cat("Begin IL2RA_T1D")
  print(Sys.time())
  IL2RA_T1D <- T1D[hm_chrom == "10"]
  IL2RA_T1D <- IL2RA_T1D[hm_pos >= 6052652 - mbp & hm_pos < 6104288 + mbp]
  saveRDS(IL2RA_T1D, file = "data/export_cis_sumstats/IL2RA_T1D.rds")

  cat("Begin IL2RB_T1D")
  print(Sys.time())
  IL2RB_T1D <- T1D[hm_chrom == "22"]
  IL2RB_T1D <- IL2RB_T1D[hm_pos >= 37521878  - mbp & hm_pos < 37571094  + mbp]
  saveRDS(IL2RB_T1D, file = "data/export_cis_sumstats/IL2RB_T1D.rds")

  cat("Begin IL6R_T1D")
  print(Sys.time())
  IL6R_T1D <- T1D[hm_chrom == "1"]
  IL6R_T1D <- IL6R_T1D[hm_pos >=  154377669 - mbp & hm_pos < 154441926 + mbp]
  saveRDS(IL6R_T1D, file = "data/export_cis_sumstats/IL6R_T1D.rds")

  cat("Begin IL6ST_T1D")
  print(Sys.time())
  IL6ST_T1D <- T1D[hm_chrom == "5"]
  IL6ST_T1D <- IL6ST_T1D[hm_pos >=  55230923 - mbp & hm_pos <  55290821 + mbp]
  saveRDS(IL6ST_T1D, file = "data/export_cis_sumstats/IL6ST_T1D.rds")

  cat("Begin TYK2_T1D")
  print(Sys.time())
  TYK2_T1D <- T1D[hm_chrom == "19"]
  TYK2_T1D <- TYK2_T1D[hm_pos >=  10461209 - mbp & hm_pos < 10491352  + mbp]
  saveRDS(TYK2_T1D, file = "data/export_cis_sumstats/TYK2_T1D.rds")

  cat("Begin JAK2_T1D")
  print(Sys.time())
  JAK2_T1D <- T1D[hm_chrom == "9"]
  JAK2_T1D <- JAK2_T1D[hm_pos >=  4985033 - mbp & hm_pos < 5128183 + mbp]
  saveRDS(JAK2_T1D, file = "data/export_cis_sumstats/JAK2_T1D.rds")

  cat("Begin JAK3_T1D")
  print(Sys.time())
  JAK3_T1D <- T1D[hm_chrom == "19"]
  JAK3_T1D <- JAK3_T1D[hm_pos >=  17935589 - mbp & hm_pos < 17958880  + mbp]
  saveRDS(JAK3_T1D, file = "data/export_cis_sumstats/JAK3_T1D.rds")


  rm(list=ls())
  gc()

  #4 read EQTL data-----------------------------------------


  #Get eQTL data from eQTLgen phase I (in unix environment)

  # All data below can be downloaded at page https://www.eqtlgen.org/cis-eqtls.html
  # Downloads -> Significant cis-eQTLs:
  # wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gzwget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
  # -> see below for decompression and renaming
  #
  # Downloads -> Full cis-eQTL summary statistics
  # wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
  # -> see below for decompression and renaming#

  # Allele frequencies
  # - > Download the allele frequencies for all tested SNPs here:
  #     Allele frequencies based on 26,609 eQTLGen samples (excluding FHS)
  # wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz
  # -> see below for decompression and renaming
  #
  #
  # #decompress and rename data
  # zcat -d 2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz > cis-EQTL-FDR0.05.txt
  # zcat -d 2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz > cis-EQTL-full.txt
  # zcat -d 2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz > cis-EQTL-AF.txt
  # move all files to data/import/

#Import data to R
  eqtl_full <- fread("data/import/cis-EQTL-full.txt")
  eqtl_AF   <- fread("data/import/cis-EQTL-AF.txt")

  #Change to table_df
  eqtl_full <- lazy_dt(eqtl_full)
  eqtl_AF  <- lazy_dt(eqtl_AF)

  #merge allele frequencies and eqtl data
  eqtl_full_AF <- left_join(eqtl_full, eqtl_AF)

  #remove original files
  rm(eqtl_AF, eqtl_full)

  #change to data table
  eqtl_full_AF <- as.data.table(eqtl_full_AF)

  #5 Filter areas of interest within the EQTL data-----------------------------------------

  #Filter only areas of interest with the genome
  eqtl_IL2RA  <- eqtl_full_AF[GeneSymbol == "IL2RA"]
  eqtl_IL2RB  <- eqtl_full_AF[GeneSymbol == "IL2RB"]
  eqtl_IL6R   <- eqtl_full_AF[GeneSymbol == "IL6R"]
  eqtl_JAK2   <- eqtl_full_AF[GeneSymbol == "JAK2"]
  eqtl_JAK3   <- eqtl_full_AF[GeneSymbol == "JAK3"]
  eqtl_TYK2   <- eqtl_full_AF[GeneSymbol == "TYK2"]
  rm(eqtl_full_AF)

  #Save output files
  saveRDS(eqtl_IL2RA,  file = "data/export_cis_sumstats/eqtl_IL2RA-Gather.rds")
  saveRDS(eqtl_IL2RB,  file = "data/export_cis_sumstats/eqtl_IL2RB-Gather.rds")
  saveRDS(eqtl_IL6R,   file = "data/export_cis_sumstats/eqtl_IL6R-Gather.rds")
  saveRDS(eqtl_JAK2,   file = "data/export_cis_sumstats/eqtl_JAK2-Gather.rds")
  saveRDS(eqtl_JAK3,   file = "data/export_cis_sumstats/eqtl_JAK3-Gather.rds")
  saveRDS(eqtl_TYK2,   file = "data/export_cis_sumstats/eqtl_TYK2-Gather.rds")


  cat("All Done with parsing")
  print(Sys.time())
