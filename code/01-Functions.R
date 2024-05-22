#1. Processing before harmonization of data------------------------------------------------------
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

check_conflicts_rsid <- function(data, effectAllele.x, effectAllele.y) {
  filter(data, rsids.x != rsids.y |
           is.na(rsids.x) & !is.na(rsids.y) |
           !is.na(rsids.x) &  is.na(rsids.y)) %>%
    mutate(rsids.x = factor(rsids.x),
           rsids.y = factor(rsids.y)) %>%
    select(rsids.x, rsids.y) %>%
    summary()
}

check_conflicts_ea_oa <- function(data, effectAllele.x, effectAllele.y) {
  filter(data, effectAllele.x != effectAllele.y |
           is.na(effectAllele.x) & !is.na(effectAllele.y) |
           !is.na(effectAllele.x) &  is.na(effectAllele.y) |
           otherAllele.x != otherAllele.y |
           is.na(otherAllele.x) & !is.na(otherAllele.y) |
           !is.na(otherAllele.x) &  is.na(otherAllele.y))
}

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

get_n_conflicts  <- function(data, AlleleB=AlleleB, AssessedAllele=AssessedAllele) {
  data %>%
    filter(AlleleB != AssessedAllele) %>%
    pull(SNP) %>% length()
}

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

#2. Colocaization analysis for single cells-----------------------------------------

# The Main function
get_coloc_hypothesis_tables <- function(exposure_gene, exposure_celltype) {

  data <- import_coloc_data(exposure_gene=exposure_gene, exposure_celltype=exposure_celltype)
  coloc_results <- analyze_coloc_celltypes(data)
  h_table <- coloc_results$coloc_results$summary
  h_table <- unname(h_table)
  h_table <- suppressMessages(h_table)
  h_table <- suppressWarnings(h_table)
  attributes(h_table) <- NULL
  h_table

}

# Supporting functions----------------------------------------------------------#
analyze_coloc_celltypes <- function(data) {

  data <- data %>%
    filter(!is.na(eaf.outcome))

  D1 <- list(
    type = "quant", # quantitative trait
    pvalues = data$pval.exposure,
    N = data$samplesize.exposure, #this is WRONG, but I'll let it be for now
    MAF = ifelse(!is.na(data$eaf.exposure),
                 ifelse(data$eaf.exposure < 0.5, data$eaf.exposure, 1  - data$eaf.exposure),
                 ifelse(data$eaf.outcome < 0.5, data$eaf.outcome, 1  - data$eaf.outcome)),
    position = data$pos.outcome,
    phenotype = data$exposure,
    snp = data$SNP,
    sdY = 1)

  D2 <- list(
    type = "cc", # case-control trait
    pvalues = data$pval.outcome,
    N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
    s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
    MAF = ifelse(data$eaf.outcome < 0.5, data$eaf.outcome, 1  - data$eaf.outcome),
    position = data$pos.outcome,
    phenotype = data$exposure,
    snp = data$SNP)

  # check_dataset_D1 <- check_dataset(D1)
  # check_dataset_D2 <- check_dataset(D2)

  coloc_D1D2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  sensitivity_coloc_D1D2 <- sensitivity_custom(coloc_D1D2, "H4 > 0.7",
                                               label1 = data$exposure[1],
                                               label2  = data$outcome[1])

  list(coloc_results = coloc_D1D2,
       coloc_sensitivity = sensitivity_coloc_D1D2)
}
import_coloc_data <- function(exposure_gene, exposure_celltype) {
  import_data(
    paste0(
      here::here("data/export_celltypes_harmonized/"),
      exposure_gene,
      "_",
      exposure_celltype,
      "_TwoSampleMR.csv"
    ),
    paste0(
      here::here("data/export_harmonization/"),
      exposure_gene,
      "_T1D_TwoSampleMR.csv"
    ),
    has_beta = FALSE
  )
}

#Custom colocalization sensitivity function to get add labels
#(tweaked by Jake Linm see email on May 14, 2024)
sensitivity_custom <- function(obj,rule="",label1="trait 1", label2="trait 2",
                               dataset1=NULL,dataset2=NULL,
                               npoints=100,doplot=TRUE,plot.manhattans=TRUE,
                               preserve.par=FALSE,
                               row=1) {
  stopifnot("list" %in% class(obj))
  stopifnot("priors" %in% names(obj))
  stopifnot("summary" %in% names(obj))
  if(rule=="")
    stop("please supply a rule to define colocalisation, eg 'H4 > thr' where thr is some probability of H4 that you accept as colocalisation")
  rule.init <- rule
  rule <- gsub("(H.)","PP.\\1.abf",rule,perl=TRUE)

  ## massage results object
  results=obj$results
  ## multiple signals?
  multiple=FALSE
  if(is.data.table(obj$summary)) { # we're not in coloc.abf anymore
    if(!(row %in% 1:nrow(obj$summary)))
      stop("row must be between 1 and ",nrow(obj$summary))
    pp <- unlist(c(obj$summary[row,grep("PP|nsnp",names(obj$summary)),with=FALSE]))
    if(paste0("SNP.PP.H4.row",row) %in% names(results)) {
      multiple=TRUE
      results[["SNP.PP.H4"]]  <- results[[paste0("SNP.PP.H4.row",row)]]
    }
    if(paste0("z.df1.row",row) %in% names(results)) { # might be passed here or in separate dataset objects
      results[["z.df1"]]  <- results[[paste0("z.df1.row",row)]]
      results[["z.df2"]]  <- results[[paste0("z.df2.row",row)]]
    } else {
      pp <- unlist(c(obj$summary[row,grep("PP|nsnp",names(obj$summary)),with=FALSE]))
    }
  } else {
    pp <- obj$summary
  }
  ## need to add z score from datasets?
  if(!is.null(dataset1) && !is.null(dataset2)) {
    df1=with(dataset1,data.table(snp=snp,position=position,z.df1=beta/sqrt(varbeta)))
    df2=with(dataset2,data.table(snp=snp,position=position,z.df2=beta/sqrt(varbeta)))
    df=merge(df1,df2,by=c("snp","position"),all=TRUE)
    results=merge(results,df,by="snp")
  }

  p12 <- obj$priors["p12"]
  p1 <- obj$priors["p1"]
  p2 <- obj$priors["p2"]
  check <- function(pp) { with(as.list(pp),eval(parse(text=rule))) }
  pass.init <- check(pp)
  message("Results ",if(check(pp)) { "pass" } else { "fail" }, " decision rule ",rule.init)

  testp12 <- 10^seq(log10(p1*p2),log10(min(p1,p1)),length.out=npoints)
  testH <- prior.snp2hyp(pp["nsnps"],p12=testp12,p1=p1,p2=p2)
  testpp <- as.data.frame(prior.adjust(summ=pp,newp12=testp12,p1=p1,p2=p2,p12=p12))
  colnames(testpp) <- gsub("(H.)","PP.\\1.abf",colnames(testpp),perl=TRUE)
  pass <- check(testpp)
  w <- which(pass)

  if(doplot) {
    H <- as.character(0:4)
    palette(c("#ffffffff",viridis(5,alpha=1)[-1]))
    op <- par('mfcol', 'mar', 'mfrow','mar','mgp','las','tck')
    on.exit(par(op))
    if(!preserve.par) {
      if(plot.manhattans)
        layout(mat = matrix(1:4,2,2),
               heights = c(1, 1), # Heights of the two rows
               widths = c(2, 3)) # Widths of the two columns
      else
        par(mfcol=c(1,2))
    }
    par(mar = c(3, 3, 2, 1) # Dist' from plot to side of page
        ,mgp = c(2, 0.4, 0) # Dist' plot to label
        ,las = 1 # Rotate y-axis text
        ,tck = -.01 # Reduce tick length
    )
    if(plot.manhattans) {
      if(!("z.df1" %in% colnames(results)) || !("z.df2" %in% colnames(results)))
        stop("please supply dataset1, dataset2, if you want to view Manhattan plots. Otherwise set plot.manhattans=FALSE.")
      manh.plot(results,1)
      title(main=paste(label1, if(multiple) { paste("row",row) } else { "" }))
      manh.plot(results,2)
      title(main=paste(label2, if(multiple) { paste("row",row) } else { "" }))
    }
    m <- list(testH,as.matrix(testpp))
    ti <-   list("Prior probabilities", "Posterior probabilities")
    for(i in 1:2) {
      ym <- if(i==1) { max(m[[i]][,-1]) } else { max(m[[i]]) }
      matplot(testp12,m[[i]],log="x",xlab="p12",ylab="Prob",
              type="b",
              bg = 1:5, # Fill colour
              pch = 21, # Shape: circles that can filed
              col="gray20",
              frame.plot = FALSE, # Remove the frame
              panel.first = abline(h = seq(0, 1, 0.2), col = "grey80"),
              ylim=c(0,ym))
      title(main=ti[[i]],adj=0)
      title(sub=paste("shaded region:",rule.init),adj=0)
      ## title(main=paste("Acceptance rule (shaded region):",rule.init))
      ## legend("topleft",pch=rep(21,5),pt.bg=1:5,legend=paste0("H",0:4))
      if(i==1)
        legend("left",inset=c(0.1,0),bg="white",pch=rep(21,5),pt.bg=1:5,pt.cex=2,legend=paste0("H",0:4))
      abline(v=p12,lty="dashed",col="gray")
      text(p12,0.5,"results",srt=90,col="gray40")
      if(any(pass))
        rect(xleft=testp12[min(w)],ybottom=0,
             xright=testp12[max(w)],ytop=1,
             col=rgb(0,1,0,alpha=0.1), border="green")
      ## add text showing rule
      ## mtext(paste("shaded region:",rule.init),side=3,adj=1)
    }
  }

  invisible(cbind(testpp,p12=testp12,pass=pass))
}

#Supporting functions not exported by coloc
prior.adjust <- function(summ,newp12,p1=1e-4,p2=1e-4,p12=1e-6) {
  if(is.list(summ) && "summary" %in% names(summ))
    summ <- summ$summary
  if(!identical(names(summ), c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")))
    stop("not a coloc summary vector")
  ## back calculate likelihoods
  f <- function(p12)
    prior.snp2hyp(summ["nsnps"],p12=p12,p1=p1,p2=p2)
  pr1 <- f(newp12)
  pr0 <- matrix(f(p12),nrow=nrow(pr1),ncol=ncol(pr1),byrow=TRUE)
  ## make everything conformable
  ## if(is.matrix(summ) && nrow(summ)==1) summ <- as.vector(summ)
  ## if(nrow(pr1)==1) pr1 <- as.vector(pr1)
  ## if(nrow(pr1)>1) pr1 <- t(pr1)
  newpp <- matrix(summ[-1],nrow=nrow(pr1),ncol=ncol(pr1),byrow=TRUE) * pr1/pr0 # prop to, not equal to
  newpp/rowSums(newpp)
}


prior.snp2hyp <- function(nsnp,p12=1e-6,p1=1e-4,p2=1e-4) {
  if(any(p12<p1*p2) || any(p12 > p1) || any(p12 > p2))
    return(NULL)
  tmp <- cbind(nsnp * p1,
               nsnp * p2,
               nsnp * (nsnp-1) * p1 * p2,
               nsnp * p12)
  tmp <- cbind(1-rowSums(tmp),tmp)
  ## if(nrow(tmp)==1) {
  ##     tmp <- c(tmp)
  ##     names(tmp) <- paste0("H",0:4)
  ## } else
  colnames(tmp) <- paste0("H",0:4)
  tmp
}

manh.plot <- function(df,wh,
                      position=if("position" %in% names(df)) {
                        df$position
                      } else {
                        1:nrow(df)
                      }) {
  znm <- if(wh==1) { "z.df1" } else {"z.df2" }
  ## print(znm)
  ## print(head(df))
  logp <- - ( pnorm(-abs(df[[znm]]),log.p=TRUE) + log(2) ) / log(10)
  ## mycol <- ifelse(A$snp %in% nCV, "red","black")
  Pal <- colorRampPalette(c('white','blue'))

  ##This adds a column of color values
  ## based on the y values
  Col <- Pal(100)[ceiling(100*df$SNP.PP.H4)]
  plot(position,logp,col="gray20",
       bg = Col, # Fill colour
       pch = 21, # Shape: circles that can filed
       frame.plot = FALSE, # Remove the frame
       xlab=if("position" %in% names(df)) {
         "Chromosome position"
       } else {
         "SNP number"
       },
       ylab="-log10(p)",
       xaxt='n')
  ## main=paste("Trait",wh))
  axis(side=1,labels=FALSE)
}


#Not dependent, but good for testing purposes--------------------------------#

get_coloc_for_celltypes_simple <- function(exposure_gene, exposure_celltype) {

  data <- import_coloc_data(exposure_gene=exposure_gene, exposure_celltype=exposure_celltype)
  res <- analyze_coloc_celltypes(data)
  return(res$coloc_results)
}
# coloc_IL2RA_cd8et <- get_coloc_for_celltypes_simple("IL2RA", "cd8et")
#end coloc------------------------------------------------------------------------#



#3. other functions--------------------------------------------------------------------
#Function to import data analysis data

import_data <- function(exposure_data_file, outcome_data_file, has_beta, min_p = 1e-400) {
  exposure_data <- read_exposure_data(exposure_data_file, sep=",", min_pval = min_p)
  outcome_data  <- read_outcome_data(outcome_data_file, sep=",", min_pval = min_p)

  if(has_beta) {
    harmonized_data <- harmonise_data(exposure_data, outcome_data)
  } else {
    harmonized_data <- harmonise_data(exposure_data, outcome_data) %>%
      rename(zscore.exposure = beta.exposure)
  }

  return(harmonized_data)
}

#Function to produce manhattan plots.
# Error: "Warning message:
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
        position = data$pos.outcome,
        snp = data$SNP,
        sdY = 1)

      D2 <- list(
        type = "cc", # case-control trait
        pvalues = data$pval.outcome,
        N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
        s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
        MAF = data$maf.outcome,  #eqtl used here in purpose
        position = data$pos.outcome,
        snp = data$SNP)

    } else if(anyNA(data$maf.exposure) == FALSE & anyNA(data$maf.outcome)) {
        D1 <- list(
          type = "quant", # quantitative trait
          pvalues = data$pval.exposure,
          N = data$samplesize.exposure,
          MAF = data$maf.exposure,
          position = data$pos.outcome,
          snp = data$SNP,
          sdY = 1)

    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf.exposure,  #eqtl used here in purpose
      position = data$pos.outcome,
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
      position = data$pos.outcome,
      snp = data$SNP,
      sdY = 1)


    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf,
      position = data$pos.outcome,
      snp = data$SNP)

    check_dataset_D1 <- check_dataset(D1)
    check_dataset_D2 <- check_dataset(D2)
    coloc_D1D2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    sensitivity_coloc_D1D2 <- sensitivity_custom(coloc_D1D2, "H4 > 0.7",
                                                 label1 = data$exposure[1],
                                                 label2  = data$outcome[1])
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
      position = data$pos.outcome,
      snp = data$SNP,
      sdY = 1)

    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      beta = IL6R$beta.outcome,
      varbeta = IL6R$se.outcome^2,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf.outcome,  #eqtl used here in purpose
      position = data$pos.outcome,
      snp = data$SNP)

  } else if(anyNA(data$maf.exposure) == FALSE & anyNA(data$maf.outcome)) {
    D1 <- list(
      type = "quant", # quantitative trait
      pvalues = data$pval.exposure,
      N = data$samplesize.exposure,
      MAF = data$maf.exposure,
      position = data$pos.outcome,
      snp = data$SNP,
      sdY = 1)

    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf.exposure,  #eqtl used here in purpose
      position = data$pos.outcome,
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
      position = data$pos.outcome,
      snp = data$SNP,
      sdY = 1)


    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf,
      position = data$pos.outcome,
      snp = data$SNP)

    check_dataset_D1 <- check_dataset(D1)
    check_dataset_D2 <- check_dataset(D2)
    coloc_D1D2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    sensitivity_coloc_D1D2 <- sensitivity(coloc_D1D2, "H4 > 0.7")
  }

  return(list("test_D1" = check_dataset_D1, "test_D2" = check_dataset_D2,
              "coloc_D1D2" = coloc_D1D2, "sensitivity" = sensitivity_coloc_D1D2))
}


analyze_coloc_pqtl_nosensitivitity <- function(data) {
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
      position = data$pos.outcome,
      snp = data$SNP,
      sdY = 1)

    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      beta = IL6R$beta.outcome,
      varbeta = IL6R$se.outcome^2,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf.outcome,  #eqtl used here in purpose
      position = data$pos.outcome,
      snp = data$SNP)

  } else if(anyNA(data$maf.exposure) == FALSE & anyNA(data$maf.outcome)) {
    D1 <- list(
      type = "quant", # quantitative trait
      pvalues = data$pval.exposure,
      N = data$samplesize.exposure,
      MAF = data$maf.exposure,
      position = data$pos.outcome,
      snp = data$SNP,
      sdY = 1)

    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf.exposure,  #eqtl used here in purpose
      position = data$pos.outcome,
      snp = data$SNP)


    check_dataset_D1 <- check_dataset(D1)
    check_dataset_D2 <- check_dataset(D2)
    coloc_D1D2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    # sensitivity_coloc_D1D2 <- sensitivity(coloc_D1D2, "H4 > 0.7")
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
      position = data$pos.outcome,
      snp = data$SNP,
      sdY = 1)


    D2 <- list(
      type = "cc", # case-control trait
      pvalues = data$pval.outcome,
      N = 18942+501638, # Case-control study (Chiou et al. 2021 Nature)
      s = 18942/(18942+501638), # N_case/(N_case+ N_ctrl)
      MAF = data$maf,
      position = data$pos.outcome,
      snp = data$SNP)

    check_dataset_D1 <- check_dataset(D1)
    check_dataset_D2 <- check_dataset(D2)
    coloc_D1D2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    # sensitivity_coloc_D1D2 <- sensitivity(coloc_D1D2, "H4 > 0.7")
  }

  return(list("test_D1" = check_dataset_D1, "test_D2" = check_dataset_D2,
              "coloc_D1D2" = coloc_D1D2))
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

  d = data

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
    geom_label_repel(aes(label = lab), min.segment.length = 0,
                    box.padding = 0.5, max.overlaps = Inf, col = "black", alpha = 0.8) +
    geom_label_repel(aes(label = lab2), min.segment.length = 0,
                    box.padding = 0.5, max.overlaps = Inf, col = "black", alpha = 0.8) +
    scale_x_continuous(expand = c(0,0),
                       limits = c(fig_start/1000, fig_end/1000)) +
    scale_colour_manual(name = bquote(italic(R)^2~.(paste0(" with ", leadSNP_trait1))),
      values = c("0.8-1.0" = "red",
                 "0.6-0.8" = "darkorange1",
                 "0.4-0.6" = "green1",
                 "0.2-0.4" = "skyblue1",
                 "0.0-0.2" = "navyblue",
                 "NA" = "grey"),
      labels = c(rev(levels(d$LD)), "Not available"),
      breaks = c(rev(levels(d$LD)), "Not available"),  # added this line to reverse the legend
      drop = FALSE) +
    theme_bw() +
    ggtitle(d$exposure[1]) +
    theme(plot.title = element_text(hjust = 0.5))

  #Figure for trait 2 (outcome)
  p2 <- ggplot(d, aes(x = pos.outcome/1000, y = -log10(pval.outcome), col = LD)) +
    geom_point() +
    labs(x = "", y = bquote(-log[10](italic(p)))) +
    geom_label_repel(aes(label = lab), min.segment.length = 0,
                    box.padding = 0.5, max.overlaps = Inf, col = "black", alpha = 0.8) +
    geom_label_repel(aes(label = lab2), min.segment.length = 0,
                    box.padding = 0.5, max.overlaps = Inf, col = "black", alpha = 0.8) +
    scale_x_continuous(expand = c(0,0),
                       limits = c(fig_start/1000, fig_end/1000)) +
    scale_colour_manual(name = bquote(italic(R)^2~.(paste0(" with ", leadSNP_trait2))),
      values = c("0.8-1.0" = "red",
                 "0.6-0.8" = "darkorange1",
                 "0.4-0.6" = "green1",
                 "0.2-0.4" = "skyblue1",
                 "0.0-0.2" = "navyblue",
                 "NA" = "grey"),
      labels = c(rev(levels(d$LD)), "Not available"),
      breaks = c(rev(levels(d$LD)), "Not available"),  # added this line to reverse the legend
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
                       hjust = 1, vjust = -1.25, size = 1.5) +
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



#This is just the same function as above, but I commented out the lead SNP for T1D
#one unnecessary label for a graph for IL6R

plot_coloc_IL6R <- function(genechr,
                       hg38_start = fig_start - gene_window,
                       hg38_end   = fig_end + gene_window,
                       fig_start,
                       fig_end,
                       gene_window = 100000,
                       coloc_results, data, D1, D2, LDmat) {

  d = data

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
    geom_label_repel(aes(label = lab), min.segment.length = 0,
                     box.padding = 0.5, max.overlaps = Inf, col = "black", alpha = 0.8) +
    # geom_label_repel(aes(label = lab2), min.segment.length = 0,
    #                  box.padding = 0.5, max.overlaps = Inf, col = "black", alpha = 0.8) +
    scale_x_continuous(expand = c(0,0),
                       limits = c(fig_start/1000, fig_end/1000)) +
    scale_colour_manual(name = bquote(italic(R)^2~.(paste0(" with ", leadSNP_trait1))),
      values = c("0.8-1.0" = "red",
                 "0.6-0.8" = "darkorange1",
                 "0.4-0.6" = "green1",
                 "0.2-0.4" = "skyblue1",
                 "0.0-0.2" = "navyblue",
                 "NA" = "grey"),
      labels = c(rev(levels(d$LD)), "Not available"),
      breaks = c(rev(levels(d$LD)), "Not available"),  # added this line to reverse the legend
      drop = FALSE) +
    theme_bw() +
    ggtitle(d$exposure[1]) +
    theme(plot.title = element_text(hjust = 0.5))

  #Figure for trait 2 (outcome)
  p2 <- ggplot(d, aes(x = pos.outcome/1000, y = -log10(pval.outcome), col = LD)) +
    geom_point() +
    labs(x = "", y = bquote(-log[10](italic(p)))) +
    geom_label_repel(aes(label = lab), min.segment.length = 0,
                     box.padding = 0.5, max.overlaps = Inf, col = "black", alpha = 0.8) +
    # geom_label_repel(aes(label = lab2), min.segment.length = 0,
    #                  box.padding = 0.5, max.overlaps = Inf, col = "black", alpha = 0.8) +
    scale_x_continuous(expand = c(0,0),
                       limits = c(fig_start/1000, fig_end/1000)) +
    scale_colour_manual(name = bquote(italic(R)^2~.(paste0(" with ", leadSNP_trait2))),
      values = c("0.8-1.0" = "red",
                 "0.6-0.8" = "darkorange1",
                 "0.4-0.6" = "green1",
                 "0.2-0.4" = "skyblue1",
                 "0.0-0.2" = "navyblue",
                 "NA" = "grey"),
      labels = c(rev(levels(d$LD)), "Not available"),
      breaks = c(rev(levels(d$LD)), "Not available"),  # added this line to reverse the legend
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
                       hjust = 1, vjust = -1.25, size = 1.5) +
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




get_ld_matrix <- function(target) {
  # Access the SNP column directly from the harmonized_dfs list
  target_variants = harmonized_dfs[[target]]$SNP

  ld_matrix(
    variants = target_variants,
    pop = "EUR",
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = plink_path,
    with_alleles = FALSE)
}

trim_matrix <- function(input_matrix, reference_data) {
  input_matrix = input_matrix[rownames(input_matrix) %in% reference_data$SNP, ]
  input_matrix =  input_matrix[, colnames(input_matrix) %in% reference_data$SNP]
  return(input_matrix)
}



plot_coloc_2 <- function(genechr,
                       hg38_start = fig_start - gene_window,
                       hg38_end   = fig_end + gene_window,
                       fig_start,
                       fig_end,
                       gene_window = 100000,
                       coloc_results, data, D1, D2, LDmat) {

  d = data

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
    geom_label_repel(aes(label = lab), min.segment.length = 0,
                     box.padding = 0.5, max.overlaps = Inf, col = "black", alpha = 0.8) +
    geom_label_repel(aes(label = lab2), min.segment.length = 0,
                     box.padding = 0.5, max.overlaps = Inf, col = "black", alpha = 0.8) +
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
    geom_label_repel(aes(label = lab), min.segment.length = 0,
                     box.padding = 0.5, max.overlaps = Inf, col = "black", alpha = 0.8) +
    geom_label_repel(aes(label = lab2), min.segment.length = 0,
                     box.padding = 0.5, max.overlaps = Inf, col = "black", alpha = 0.8) +
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
    geom_text_repel(data = plot_lab,
                       aes(x = value/1000, y = y_value, label = GENE_NAME),
                       size = 1.5) +
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


sig2 <- function(x) signif(x, digits = 2)



