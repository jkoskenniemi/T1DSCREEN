
#Function to import data analysis data

import_data <- function(exposure_data_file, outcome_data_file, has_beta) {
  exposure_data <- read_exposure_data(exposure_data_file, sep=",")
  outcome_data  <- read_outcome_data(outcome_data_file, sep=",")

  if(has_beta) {
    harmonized_data <- harmonise_data(exposure_data, outcome_data)
  } else {
    harmonized_data <- harmonise_data(exposure_data, outcome_data) %>%
      rename(zscore.exposure = beta.exposure)
  }

  return(harmonized_data)
}

#Function to produce manhattan plots.
# Error
# "Warning message:
# Removed 1 rows containing missing values (`geom_text()`). "
# usually means that you have specified incorrect start and end
# data for the gene (or at least they are not compatible with
# the data that is used)
manhattan_plots <- function(data,
                            gene_start, gene_end, title) {

  manhattan_exposure <- data %>%
    ggplot(aes(x = pos.outcome/1000, y=-log10(pval.exposure))) +
    geom_point() +
    geom_hline(yintercept=8, linetype="dashed") +
    ggtitle(paste0("GWAS of ", data$exposure[1])) +
    ylab("-log10(p)") +
    xlab(NULL) +
    theme(axis.text.x = element_blank())

  manhattan_outcome <- data %>%
    ggplot() +
    geom_point(aes(x = pos.outcome/1000, y=-log10(pval.outcome))) +
    geom_hline(yintercept=8, linetype="dashed") +
    ggtitle(paste0("GWAS of ", data$outcome[1])) +
    ylab("-log10(p)") +
    xlab("position(kbp)")

  gene_figure <- data %>%
    ggplot() +
    geom_blank() +
    geom_segment(x=gene_start/1000, xend=gene_end/1000, y=1, yend=1, linewidth = 2) +
    annotate("text", label = paste0(title),
             x= (gene_start + gene_end)/ 1000 / 2, y=1, hjust = 0.5, vjust = -0.5) +
    xlim(layer_scales(manhattan_exposure)$x$range$range) +
    ylim(0.75, 2) +
    ylab(NULL) + xlab(NULL) +
    theme(axis.text = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())

  manhattan_plot <- ggarrange(gene_figure, manhattan_exposure, manhattan_outcome,
                              heights = c(1.2, 3, 3), nrow = 3,
                              ncol = 1, align = "hv")
  return(manhattan_plot)
}


analyze_coloc_eqtl <- function(data) {
  data <- data %>%
    mutate(eaf.exposure = as.numeric(eaf.exposure)) %>%
    mutate(eaf.outcome  = as.numeric(eaf.outcome)) %>%
    mutate(maf.exposure  =
             ifelse(eaf.exposure  < 0.5, eaf.exposure,
                    1 - eaf.exposure)) %>%
    mutate(maf.outcome  =
             ifelse(eaf.outcome  < 0.5, eaf.outcome,
                    1 - eaf.outcome))

  data <- data %>%
    filter(!is.na(pos.outcome))

    if(anyNA(data$maf.exposure) & anyNA(data$maf.outcome) == FALSE) {
      D1 <- list(
        type = "quant", # quantitative trait
        pvalues = data$pval.exposure,
        N = data$samplesize.exposure,
        MAF = data$maf.outcome,
        pos = data$pos.outcome,
        snp = data$pos.outcome,
        sdY = 1)

      D2 <- list(
        type = "cc", # case-control trait
        pvalues = data$pval.outcome,
        N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
        s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
        MAF = data$maf.outcome,  #eqtl used here in purpose
        pos = data$pos.outcome,
        snp = data$pos.outcome)

    } else if(anyNA(data$maf.exposure) == FALSE & anyNA(data$maf.outcome)) {
        D1 <- list(
          type = "quant", # quantitative trait
          pvalues = data$pval.exposure,
          N = data$samplesize.exposure,
          MAF = data$maf.exposure,
          pos = data$pos.outcome,
          snp = data$pos.outcome,
          sdY = 1)

    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf.exposure,  #eqtl used here in purpose
      pos = data$pos.outcome,
      snp = data$pos.outcome)


    check_dataset_D1 <- check_dataset(D1)
    check_dataset_D2 <- check_dataset(D2)
    coloc_D1D2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    sensitivity_coloc_D1D2 <- sensitivity(coloc_D1D2, "H4 > 0.7")
  } else {

    data <- data %>%
      mutate(more_non_na_maf = ifelse(length(!is.na(maf.exposure)) >
                                 length(!is.na(maf.outcome)),
                               "exposure", "outcome")) %>%
      mutate(maf = ifelse(more_non_na_maf == "exposure", maf.exposure, maf.outcome)) %>%
      filter(!is.na(maf))

    D1 <- list(
      type = "quant", # quantitative trait
      pvalues = data$pval.exposure,
      N = data$samplesize.exposure,
      MAF = data$maf,
      pos = data$pos.outcome,
      snp = data$pos.outcome,
      sdY = 1)


    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf,
      pos = data$pos.outcome,
      snp = data$pos.outcome)

    check_dataset_D1 <- check_dataset(D1)
    check_dataset_D2 <- check_dataset(D2)
    coloc_D1D2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    sensitivity_coloc_D1D2 <- sensitivity(coloc_D1D2, "H4 > 0.7")
  }

  return(list("test_D1" = check_dataset_D1, "test_D2" = check_dataset_D2,
              coloc = "coloc_D1D2", "sensitivity" = sensitivity_coloc_D1D2))
}



analyze_coloc_pqtl <- function(data) {
  data <- data %>%
    mutate(eaf.exposure = as.numeric(eaf.exposure)) %>%
    mutate(eaf.outcome  = as.numeric(eaf.outcome)) %>%
    mutate(maf.exposure  =
             ifelse(eaf.exposure  < 0.5, eaf.exposure,
                    1 - eaf.exposure)) %>%
    mutate(maf.outcome  =
             ifelse(eaf.outcome  < 0.5, eaf.outcome,
                    1 - eaf.outcome))

  data <- data %>%
    filter(!is.na(pos.outcome))

  if(anyNA(data$maf.exposure) & anyNA(data$maf.outcome) == FALSE) {
    D1 <- list(
      type = "quant", # quantitative trait
      pvalues = data$pval.exposure,
      beta = data$beta.exposure,
      varBeta = data$beta.exposure^2,
      N = data$samplesize.exposure,
      MAF = data$maf.outcome,
      pos = data$pos.outcome,
      snp = data$pos.outcome,
      sdY = 1)

    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      beta = IL6R$beta.outcome,
      varbeta = IL6R$se.outcome^2,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf.outcome,  #eqtl used here in purpose
      pos = data$pos.outcome,
      snp = data$pos.outcome)

  } else if(anyNA(data$maf.exposure) == FALSE & anyNA(data$maf.outcome)) {
    D1 <- list(
      type = "quant", # quantitative trait
      pvalues = data$pval.exposure,
      N = data$samplesize.exposure,
      MAF = data$maf.exposure,
      pos = data$pos.outcome,
      snp = data$pos.outcome,
      sdY = 1)

    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf.exposure,  #eqtl used here in purpose
      pos = data$pos.outcome,
      snp = data$pos.outcome)


    check_dataset_D1 <- check_dataset(D1)
    check_dataset_D2 <- check_dataset(D2)
    coloc_D1D2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    sensitivity_coloc_D1D2 <- sensitivity(coloc_D1D2, "H4 > 0.7")
  } else {

    data <- data %>%
      mutate(more_non_na_maf = ifelse(length(!is.na(maf.exposure)) >
                                        length(!is.na(maf.outcome)),
                                      "exposure", "outcome")) %>%
      mutate(maf = ifelse(more_non_na_maf == "exposure", maf.exposure, maf.outcome)) %>%
      filter(!is.na(maf))

    D1 <- list(
      type = "quant", # quantitative trait
      pvalues = data$pval.exposure,
      N = data$samplesize.exposure,
      MAF = data$maf,
      pos = data$pos.outcome,
      snp = data$pos.outcome,
      sdY = 1)


    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf,
      pos = data$pos.outcome,
      snp = data$pos.outcome)

    check_dataset_D1 <- check_dataset(D1)
    check_dataset_D2 <- check_dataset(D2)
    coloc_D1D2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    sensitivity_coloc_D1D2 <- sensitivity(coloc_D1D2, "H4 > 0.7")
  }

  return(list("test_D1" = check_dataset_D1, "test_D2" = check_dataset_D2,
              coloc = "coloc_D1D2", "sensitivity" = sensitivity_coloc_D1D2))
}


analyze_coloc_eqtl_pull_top_snp <- function(data) {
  data <- data %>%
    mutate(eaf.exposure = as.numeric(eaf.exposure)) %>%
    mutate(eaf.outcome  = as.numeric(eaf.outcome)) %>%
    mutate(maf.exposure  =
             ifelse(eaf.exposure  < 0.5, eaf.exposure,
                    1 - eaf.exposure)) %>%
    mutate(maf.outcome  =
             ifelse(eaf.outcome  < 0.5, eaf.outcome,
                    1 - eaf.outcome))

  data <- data %>%
    filter(!is.na(pos.outcome))

  if(anyNA(data$maf.exposure) & anyNA(data$maf.outcome) == FALSE) {
    D1 <- list(
      type = "quant", # quantitative trait
      pvalues = data$pval.exposure,
      N = data$samplesize.exposure,
      MAF = data$maf.outcome,
      pos = data$pos.outcome,
      snp = data$pos.outcome,
      sdY = 1)

    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf.outcome,  #eqtl used here in purpose
      pos = data$pos.outcome,
      snp = data$pos.outcome)

  } else if(anyNA(data$maf.exposure) == FALSE & anyNA(data$maf.outcome)) {
    D1 <- list(
      type = "quant", # quantitative trait
      pvalues = data$pval.exposure,
      N = data$samplesize.exposure,
      MAF = data$maf.exposure,
      pos = data$pos.outcome,
      snp = data$pos.outcome,
      sdY = 1)

    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf.exposure,  #eqtl used here in purpose
      pos = data$pos.outcome,
      snp = data$pos.outcome)


    coloc_D1D2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

  } else {

    data <- data %>%
      mutate(more_non_na_maf = ifelse(length(!is.na(maf.exposure)) >
                                        length(!is.na(maf.outcome)),
                                      "exposure", "outcome")) %>%
      mutate(maf = ifelse(more_non_na_maf == "exposure", maf.exposure, maf.outcome)) %>%
      filter(!is.na(maf))

    D1 <- list(
      type = "quant", # quantitative trait
      pvalues = data$pval.exposure,
      N = data$samplesize.exposure,
      MAF = data$maf,
      pos = data$pos.outcome,
      snp = data$pos.outcome,
      sdY = 1)


    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf,
      pos = data$pos.outcome,
      snp = data$pos.outcome)

    coloc_D1D2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

  }

  return(coloc_results$results[which(coloc_results$results$SNP.PP.H4==max(coloc_results$results$SNP.PP.H4)),"snp"])
}

