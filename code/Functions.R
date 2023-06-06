
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
        snp = data$SNP,
        sdY = 1)

      D2 <- list(
        type = "cc", # case-control trait
        pvalues = data$pval.outcome,
        N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
        s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
        MAF = data$maf.outcome,  #eqtl used here in purpose
        pos = data$pos.outcome,
        snp = data$SNP)

    } else if(anyNA(data$maf.exposure) == FALSE & anyNA(data$maf.outcome)) {
        D1 <- list(
          type = "quant", # quantitative trait
          pvalues = data$pval.exposure,
          N = data$samplesize.exposure,
          MAF = data$maf.exposure,
          pos = data$pos.outcome,
          snp = data$SNP,
          sdY = 1)

    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf.exposure,  #eqtl used here in purpose
      pos = data$pos.outcome,
      snp = data$SNP)


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
      snp = data$SNP,
      sdY = 1)


    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf,
      pos = data$pos.outcome,
      snp = data$SNP)

    check_dataset_D1 <- check_dataset(D1)
    check_dataset_D2 <- check_dataset(D2)
    coloc_D1D2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    sensitivity_coloc_D1D2 <- sensitivity(coloc_D1D2, "H4 > 0.7")
  }

  return(list("test_D1" = check_dataset_D1,
              "test_D2" = check_dataset_D2,
              "coloc_D1D2" = coloc_D1D2,
              "D1" = D1,
              "D2" = D2,
              "sensitivity" = sensitivity_coloc_D1D2))
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
      snp = data$SNP,
      sdY = 1)


    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf,
      pos = data$pos.outcome,
      snp = data$SNP)

    coloc_D1D2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

  }

  return(coloc_results$results[which(coloc_results$results$SNP.PP.H4==max(coloc_results$results$SNP.PP.H4)),"snp"])
}

# rm(list=ls())


#This function is based on

#This is very much a hacked version of the plots in RACER package (https://github.com/oliviasabik/RACER).
#In particular, the genetic annotation in the bottom panel is based on the "hg19" data object included in
#the RACER package.

plot_coloc <- function(genechr,
                       hg38_start = fig_start - gene_window,
                       hg38_end   = fig_end + gene_window,
                       fig_start,
                       fig_end,
                       gene_window = 100000,
                       coloc_results, data, D1, D2, LDmat) {

  coloc_results = coloc_results_IL2RA_eqtl$coloc_D1D2
  d = data
  D1 = coloc_results$D1
  D2 = coloc_results$D2
  LDmat = LDmat_IL2RA
  genome_version = "hg38"


  stopifnot(class(coloc_results)[1] == "coloc_abf")
  pacman::p_load("coloc", "data.table", "ggplot2", "ggpubr", "ggrepel","reshape2", "RACER")

  #Find the lead SNP rsid from coloc_results for trait 1
  leadSNP_trait1 <- coloc_results[["results"]][
    which.max(coloc_results[["results"]]$lABF.df1),"snp"]

  #Find the lead SNP rsid from coloc_results for trait 2
  leadSNP_trait2 <- coloc_results[["results"]][
    which.max(coloc_results[["results"]]$lABF.df2),"snp"]

  #Make a data.frame for the correlation matrix for both SNPs
  d_ld <- data.frame(rsid = colnames(LDmat),
                     LDRsq = LDmat[leadSNP_trait1,]^2)
  d_ld2 <- data.frame(rsid = colnames(LDmat), LDRsq2 = LDmat[leadSNP_trait2,]^2)

  #change from data.frame to data.table
  setDT(d_ld)
  setDT(d_ld2)
  setDT(d)


  #Color SNPs according to the R-squared with the lead SNP
  d_ld[,"LD" := cut(LDRsq, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                    labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"))]
  d_ld2[,"LD2" := cut(LDRsq2, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                      labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"))]

  #Merge LD data and data
  d <- d[d_ld, on = c("SNP" = "rsid"), nomatch = NULL]
  d <- d[d_ld2, on = c("SNP" = "rsid"), nomatch = NULL]

  #Annotate lead SNPs in both traits
  d$lab <- ifelse(d$SNP %in% leadSNP_trait1, leadSNP_trait1, "")
  d$lab2 <- ifelse(d$SNP %in% leadSNP_trait2, leadSNP_trait2, "")

  #Figure for trait 1 (exposure)
  p1 <- ggplot(d, aes(x = pos.outcome/1000, y = -log10(pval.exposure), col = LD)) +
    geom_point() +
    labs(x = "", y = bquote(-log[10](italic(p)))) +
    geom_text_repel(aes(label = lab), min.segment.length = 0,
                    box.padding = 0.5, max.overlaps = Inf, col = "black") +
    geom_text_repel(aes(label = lab2), min.segment.length = 0,
                    box.padding = 0.5, max.overlaps = Inf, col = "black") +
    scale_x_continuous(expand = c(0,0),
                       limits = c(fig_start/1000, fig_end/1000)) +
    scale_colour_manual(name = bquote(italic(R)^2~.(paste0(" with ", leadSNP_trait1))),
                        values = c("0.8-1.0" = "red",
                                   "0.6-0.8" = "darkorange1",
                                   "0.4-0.6" = "green1",
                                   "0.2-0.4" = "skyblue1",
                                   "0.0-0.2" = "navyblue",
                                   "NA" = "grey"),
                        labels = c(levels(d$LD), "Not available"),
                        drop = FALSE) +
    theme_bw() +
    ggtitle(d$exposure[1]) +
    theme(plot.title = element_text(hjust = 0.5))

  #Figure for trait 2 (outcome)
  p2 <- ggplot(d, aes(x = pos.outcome/1000, y = -log10(pval.outcome), col = LD)) +
    geom_point() +
    labs(x = "", y = bquote(-log[10](italic(p)))) +
    geom_text_repel(aes(label = lab), min.segment.length = 0,
                    box.padding = 0.5, max.overlaps = Inf, col = "black") +
    geom_text_repel(aes(label = lab2), min.segment.length = 0,
                    box.padding = 0.5, max.overlaps = Inf, col = "black") +
    scale_x_continuous(expand = c(0,0),
                       limits = c(fig_start/1000, fig_end/1000)) +
    scale_colour_manual(name = bquote(italic(R)^2~.(paste0(" with ", leadSNP_trait2))),
                        values = c("0.8-1.0" = "red",
                                   "0.6-0.8" = "darkorange1",
                                   "0.4-0.6" = "green1",
                                   "0.2-0.4" = "skyblue1",
                                   "0.0-0.2" = "navyblue",
                                   "NA" = "grey"),
                        labels = c(levels(d$LD), "Not available"),
                        drop = FALSE) +
    ggtitle(d$outcome[1]) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  # Figure 3: genes
  # add gene info
  utils::data(hg38)
  colnames(hg38) = c("GENE_ID", "CHR", "TRX_START",
                     "TRX_END", "LENGTH", "GENE_NAME",
                     "TYPE")
  gene_sub = hg38[hg38$CHR == genechr, ]
  gene_sub = gene_sub[gene_sub$TYPE == "protein_coding",]
  gene_sub = subset(gene_sub, gene_sub$TRX_START > ({{hg38_start}}))
  gene_sub = subset(gene_sub, gene_sub$TRX_END < ({{hg38_end}}))
  gene_sub = gene_sub[!duplicated(gene_sub$GENE_ID), ]
  gene_sub = gene_sub[, c(3, 4, 6)]


  #remove genes that don't stretch to the area of the figure
  gene_sub = subset(gene_sub, TRX_END > fig_start | TRX_START < fig_end)


  #If the gene starts or ends outside figure, set gene to remain in figure
  gene_sub$TRX_END   = ifelse(gene_sub$TRX_END   > fig_end,   fig_end-1,   gene_sub$TRX_END)
  gene_sub$TRX_START = ifelse(gene_sub$TRX_START < fig_start, fig_start+1, gene_sub$TRX_START)

  #Change from wide format (with both end and start per row) to long format
  #(end and start as separate observations)
  gene_sub = reshape2::melt(gene_sub, id.vars = "GENE_NAME")

  #Refresh y-axis to exclude empty rows
  gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))

  #Calculate x-axis midpoints for each gene
  plot_lab = gene_sub %>%
    group_by(GENE_NAME) %>%
    summarise(value = mean(value), y_value = mean(y_value))

  #Draw gene figure
  p_genes <- ggplot2::ggplot(gene_sub,
                             ggplot2::aes(x = value/1000, y = y_value)) +
    ggplot2::geom_line(aes(group = GENE_NAME), size = 2) +
    ggplot2::theme_bw() +
    ggplot2::geom_text(data = plot_lab,
                       aes(x = value/1000, y = y_value, label = GENE_NAME),
                       hjust = 0, vjust = -1.25, size = 1.5) +
    ggplot2::scale_x_continuous(expand = c(0,0),
                                limits = c(fig_start/1000, fig_end/1000)) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_text(color = "white", size = 28),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::xlab(paste0("Chromosome ", genechr, " Position (kB)")) +
    ggplot2::ylim(0, (max(gene_sub$y_value) + 1))

  #Combine the three figures into 1
  ggpubr::ggarrange(p1, p2, p_genes,
                    heights = c(1.5, 1.5, 1), nrow = 3,
                    ncol = 1, common.legend = TRUE, legend = "right", align = "hv")
}

# plot_coloc(genechr = 10,
#            fig_start = 5900000,
#            fig_end = 6250000,
#            gene_window = 50000,
#            coloc_results = coloc_results_IL2RA_eqtl$coloc_D1D2,
#            data = IL2RA_eqtl,
#            D1 = coloc_results_IL2RA_eqtl$D1,
#            D2 = coloc_results_IL2RA_eqtl$D2,
#            LDmat = LDmat_IL2RA)



