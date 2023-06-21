
#Load packages

#If run for the first time
# install.packages("devtools")
# library(devtools)
# install_github("phenoscanner/phenoscanner")
# library(phenoscanner)

#Load packages
pacman::p_load(tidyverse, phenoscanner)

names_gwas <- c("IFNAR2", "IL2RA", "IL2RB",
                "IL12B", "IL2RG", "IL6R",
                "IL6ST", "JAK1", "JAK2",
                "JAK3", "TYK2")

functional_variants <- map(names_gwas, ~phenoscanner(genequery = .x))
names(functional_variants) <- names_gwas

save(functional_variants, file = "data/export/phenoscanner-functional-variants.RDS")
load("data/export/phenoscanner-functional-variants.RDS")


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


save(filtered_missense, file = "data/export/phenoscanner-functional-variants-missense.csv")


save(functional_variants, file = "data/export/phenoscanner-functional-variants-rsids.RDS")



