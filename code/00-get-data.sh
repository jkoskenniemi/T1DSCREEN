#!/bin/bash

# Script to get data from DECODE, GWAS catalog and eQTLgen

# Open terminal and move to folder T1DSCREEN/data
cd T1DSCREEN/data

# Get T1D data from GWAS catalog (Chiou et al.)
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014023/harmonised/34012112-GCST90014023-EFO_0001359.h.tsv.gz
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014023/harmonised/md5sum.txt
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014023/README.txt

# Get T1D data from GWAS catalog (Robertson et al.)
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013445/GCST90013445_buildGRCh38.tsv
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013445/README.txt
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013445/md5sum.txt

# Get eQTLgen data
wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz
wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/README_cis

# Instructions for DECODE data
echo "To get DECODE data, please visit https://www.decode.com/summarydata/ and download the following files:"
echo "1. Summary data (list): https://download.decode.is/form/folder/proteomics"
echo "2. Extra annotation: https://download.decode.is/form/2021/assocvariants.annotated.txt.gz"
echo "3. Readme file: https://download.decode.is/form/2021/proteomics_readme.txt"
echo "Save the files to the data/import/ directory."

# Decompress data
zcat 34012112-GCST90014023-EFO_0001359.h.tsv.gz > 34012112-GCST90014023-EFO_0001359.h.tsv
zcat 2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz > cis-EQTL-full.txt
zcat 2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz > cis-EQTL-AF.txt

# Move data to import
mkdir -p import
mv *.tsv import/
mv *.txt import/
mv README_cis import/

# Remove all gz files
rm *.gz

# Prepare for analyses
mkdir -p export_cis_sumstats
mkdir -p export_manuscript
mkdir -p export_celltypes_harmonized
mkdir -p export_LDmat

echo "Script completed successfully."

