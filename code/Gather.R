library(data.table)
mbp <- 1000000
setwd("/scratch/project_2006579/")

# Changes Apr 25, 2023
# 1. An attempt is made to make the entire analysis workflow traceable. 
# 2. Export path is changed from /data/sumstats/ to /data/cis_sumstats/
# 3. Uncommented and modified IL2RG (protein, annotation & T1D sumstats)

#1 Read circulating protein levels---------------------------------------------

#Read annotations of Ferkingstad et al. 2021 data. This will be merged with each
#separately with each protein data sets (such as IL6R).
anno  <- fread("data/sumstats/assocvariants.annotated.txt")

#IL2RA

#Read data
cat("Start parsing IL2RA")
print(Sys.time())
IL2RA <- fread("data/sumstats/3151_6_IL2RA_IL_2_sRa.txt")
#Filter only those variants in the same chromosome 10
IL2RA <- IL2RA[Chrom == "chr10"]
#filter based on start and end Position of the gene
IL2RA <- IL2RA[Pos >= 6052652 - mbp & Pos < 6104288 + mbp]
#Filter annotation data that includes only the same chromosomes
IL2RA_anno <- anno[Chrom == "chr10"]
IL2RA_anno <- IL2RA_anno[Pos >= 6052652 - mbp & Pos < 6104288 + mbp]
saveRDS(IL2RA, file = "data/cis_sumstats/IL2RA_prot.rds")
saveRDS(IL2RA_anno, file = "data/cis_sumstats/IL2RA_anno.rds")

cat("Done parsing IL2RA, begin IL6R")
print(Sys.time())
#IL6R (similar to IL2RA)
IL6R <- fread("data/sumstats/15602_43_IL6R_IL_6_sRa.txt")
IL6R <- IL6R[Chrom == "chr1"]
IL6R <- IL6R[Pos >=  154377669 - mbp & Pos < 154441926 + mbp]
IL6R_anno <- anno[Chrom == "chr1"]
IL6R_anno <- IL6R_anno[Pos >=  154377669 - mbp & Pos < 154441926 + mbp]
saveRDS(IL6R, file = "data/cis_sumstats/IL6R_prot.rds")
saveRDS(IL6R_anno, file = "data/cis_sumstats/IL6R_anno.rds")

#IL6ST
IL6ST <- fread("data/sumstats/2620_4_IL6ST_gp130__soluble.txt")
IL6ST <- IL6ST[Chrom == "chr5"]
IL6ST <- IL6ST[Pos >=  55230923 - mbp & Pos < 55290821 + mbp]
IL6ST_anno <- anno[Chrom == "chr5"]
IL6ST_anno <- IL6ST_anno[Pos >=  55230923 - mbp & Pos <  55290821 + mbp]
saveRDS(IL6ST, file = "data/cis_sumstats/IL6ST_prot.rds")
saveRDS(IL6ST_anno, file = "data/cis_sumstats/IL6ST_anno.rds")

#TYK2
TYK2  <- fread("data/sumstats/5260_80_TYK2_TYK2.txt")
TYK2  <- TYK2 [Chrom == "chr19"]
TYK2  <- TYK2 [Pos >=  10461209 - mbp & Pos < 10491352 + mbp]
TYK2_anno <- anno[Chrom == "chr19"]
TYK2_anno <- TYK2_anno[Pos >=  10461209 - mbp & Pos < 10491352  + mbp]
saveRDS(TYK2, file = "data/cis_sumstats/TYK2_prot.rds")
saveRDS(TYK2_anno, file = "data/cis_sumstats/TYK2_anno.rds")

#JAK2
JAK2 <- fread("data/sumstats/11816_84_JAK2_JAK2.txt")
JAK2 <- JAK2[Chrom == "chr9"]
JAK2 <- JAK2[Pos >=  4985033 - mbp & Pos < 5128183 + mbp]
JAK2_anno <- anno[Chrom == "chr9"]
JAK2_anno <- JAK2_anno[Pos >=  4985033 - mbp & Pos < 5128183 + mbp]
saveRDS(JAK2, file = "data/cis_sumstats/JAK2_prot.rds")
saveRDS(JAK2_anno, file = "data/cis_sumstats/JAK2_anno.rds")

#IL12B
IL12B <- fread("data/sumstats/13733_5_IL12B_IL_12_p40.txt")
IL12B <- IL12B[Chrom == "chr5"]
IL12B <- IL12B[Pos >= 158741791  - mbp & Pos < 158757895 + mbp]
IL12B_anno <- anno[Chrom == "chr5"]
IL12B_anno <- IL12B_anno[Pos >=  158741791 - mbp & Pos < 158757895 + mbp]
saveRDS(IL12B, file = "data/cis_sumstats/IL12B_prot.rds")
saveRDS(IL12B_anno, file = "data/cis_sumstats/IL12B_anno.rds")

#IL2RG: If I remember correctly, this one was a little difficult because
#we need to filter X-chromosome. I left it out for safety at first.
IL2RG <- fread("data/sumstats/2634_2_IL2RG_IL_2_sRg.txt")
IL2RG <- IL2RG[Chrom == "chrX"]
IL2RG <- IL2RG[Pos >= 71107404  - mbp & Pos < 71111577 + mbp]
IL2RG_anno <- anno[Chrom == "chrX"]
IL2RG_anno <- IL2RG_anno[Pos >= 71107404  - mbp & Pos < 71111577 + mbp]
saveRDS(IL2RG, file = "data/cis_sumstats/IL2RB_prot.rds")
saveRDS(IL2RG_anno, file = "data/cis_sumstats/IL2RB_anno.rds")


#IL2RB
IL2RB <- fread("data/sumstats/9343_16_IL2RB_IL_2_sRb.txt")
IL2RB <- IL2RB[Chrom == "chr22"]
IL2RB <- IL2RB[Pos >= 37521878  - mbp & Pos < 37571094 + mbp]
IL2RB_anno <- anno[Chrom == "chr22"]
IL2RB_anno <- IL2RB_anno[Pos >= 37521878  - mbp & Pos < 37571094  + mbp]
saveRDS(IL2RB, file = "data/cis_sumstats/IL2RB_prot.rds")
saveRDS(IL2RB_anno, file = "data/cis_sumstats/IL2RB_anno.rds")


#CXCL10
#4141_79_CXCL10_IP_10
CXCL10 <- fread("data/sumstats/4141_79_CXCL10_IP_10.txt")
CXCL10 <- CXCL10[Chrom == "chr4"]
CXCL10 <- CXCL10[Pos >= 76021118  - mbp & Pos < 76023497 + mbp]
CXCL10_anno <- anno[Chrom == "chr4"]
CXCL10_anno <- CXCL10_anno[Pos >= 76023497  - mbp & Pos < 76023497  + mbp]
saveRDS(CXCL10, file = "data/cis_sumstats/CXCL10_prot.rds")
saveRDS(CXCL10_anno, file = "data/cis_sumstats/CXCL10_anno.rds")

#TNF alpha
#5936_53_TNF_TNF_a.txt
TNF <- fread("data/sumstats/5936_53_TNF_TNF_a.txt")
TNF <- TNF[Chrom == "chr6"]
TNF <- TNF[Pos >= 31575565  - mbp & Pos < 31578336 + mbp]

TNF_anno <- anno[Chrom == "chr6"]
TNF_anno <- TNF_anno[Pos >= 31575565  - mbp & Pos < 31578336  + mbp]
saveRDS(TNF, file = "data/cis_sumstats/TNF_prot.rds")
saveRDS(TNF_anno, file = "data/cis_sumstats/TNF_anno.rds")

#2 Read GWAS of CRP levels-----------------------------------------------------
IL6R_CRP <- fread("data/sumstats/35459240-GCST90029070-EFO_0004458.h.tsv")
IL6R_CRP <- IL6R_CRP[hm_chrom == "1"]
IL6R_CRP <- IL6R_CRP[hm_pos > 154377669 - mbp & hm_pos < 154441926 + mbp]
saveRDS(IL6R_CRP, file = "data/cis_sumstats/IL6R_CRP.rds")

cat("Done parsing proteins, begin T1D GWAS selection of protein regions")
print(Sys.time())

rm(list=ls())

#3 Read data of type 1 diabetes------------------------------------------------
T1D <- fread("data/sumstats/34012112-GCST90014023-EFO_0001359.h.tsv")

#4 Filter areas of interest in GWAS of T1D for each protein--------------------
cat("Begin IL2RA_T1D")
print(Sys.time())
IL2RA_T1D <- T1D[hm_chrom == "10"]
IL2RA_T1D <- IL2RA_T1D[hm_pos >= 6052652 - mbp & hm_pos < 6104288 + mbp]
saveRDS(IL2RA_T1D, file = "data/cis_sumstats/IL2RA_T1D.rds")

cat("Begin IL6R_T1D")
print(Sys.time())
IL6R_T1D <- T1D[hm_chrom == "1"]
IL6R_T1D <- IL6R_T1D[hm_pos >=  154377669 - mbp & hm_pos < 154441926 + mbp]
saveRDS(IL6R_T1D, file = "data/cis_sumstats/IL6R_T1D.rds")

cat("Begin IL6ST_T1D")
print(Sys.time())
IL6ST_T1D <- T1D[hm_chrom == "5"]
IL6ST_T1D <- IL6ST_T1D[hm_pos >=  55230923 - mbp & hm_pos <  55290821 + mbp]
saveRDS(IL6ST_T1D, file = "data/cis_sumstats/IL6ST_T1D.rds")

cat("Begin TYK2_T1D")
print(Sys.time())
TYK2_T1D <- T1D[hm_chrom == "19"]
TYK2_T1D <- TYK2_T1D[hm_pos >=  10461209 - mbp & hm_pos < 10491352  + mbp]
saveRDS(TYK2_T1D, file = "data/cis_sumstats/TYK2_T1D.rds")

cat("Begin JAK2_T1D")
print(Sys.time())
JAK2_T1D <- T1D[hm_chrom == "9"]
JAK2_T1D <- JAK2_T1D[hm_pos >=  4985033 - mbp & hm_pos < 5128183 + mbp]
saveRDS(JAK2_T1D, file = "data/cis_sumstats/JAK2_T1D.rds")

cat("Begin IL12B_T1D")
print(Sys.time())
IL12B_T1D <- T1D[hm_chrom == "5"]
IL12B_T1D <- IL12B_T1D[hm_pos >=  158741791 - mbp & hm_pos < 158757895 + mbp]
saveRDS(IL12B_T1D, file = "data/cis_sumstats/IL12B_T1D.rds")

cat("Begin IL2RG_T1D")
print(Sys.time())
IL2RG_T1D <- T1D[hm_chrom == "X"]
IL2RG_T1D <- IL2RG_T1D[hm_pos >= 71107404 - mbp & hm_pos < 71111577 + mbp]
saveRDS(IL2RG_T1D, file = "data/cis_sumstats/IL2RG_T1D.rds")

cat("Begin IL2RB_T1D")
print(Sys.time())
IL2RB_T1D <- T1D[hm_chrom == "22"]
IL2RB_T1D <- IL2RB_T1D[hm_pos >= 37521878  - mbp & hm_pos < 37571094  + mbp]
saveRDS(IL2RB_T1D, file = "data/cis_sumstats/IL2RB_T1D.rds")

cat("Begin CXCL10_T1D")
print(Sys.time())
CXCL10_T1D <- T1D[hm_chrom == "4"]
CXCL10_T1D <- CXCL10_T1D[hm_pos >=  76021118 - mbp & hm_pos < 76023497 + mbp]
saveRDS(CXCL10_T1D, file = "data/cis_sumstats/CXCL10_T1D.rds")

cat("Begin TNF_T1d")
print(Sys.time())
TNF_T1D <- T1D[hm_chrom == "6"]
TNF_T1D <- TNF_T1D[hm_pos >= 31575565 - mbp & hm_pos < 31578336 + mbp]
saveRDS(TNF_T1D, file = "data/cis_sumstats/TNF_T1D.rds")

cat("All Done with parsing")
print(Sys.time())
