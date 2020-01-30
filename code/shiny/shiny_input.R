
library(DT)
library(RColorBrewer)
library(strex)
library(TwoSampleMR)
library(tidyverse)

#deployApp(appDir = "code/shiny", appName = "proxymr")

load('data.RData')
#source("../functions.R")

MR_method_list <- c("mr_wald_ratio","mr_ivw","mr_ivw_fe","mr_egger_regression","mr_weighted_median")

household_GWAS_result <- c("gwas_age_bins", "gwas_tt_bins")
sex <- c("male","female")
grouping_var <- c("age_even_bins", "time_together_even_bins")

traits$shiny_description <- paste0(traits$SGG_PHESANT_ID, ": ", traits$description, " (", traits$variable_type, ")")
trait_list <- traits$shiny_description
phenotype_list <- c("raw", "phesant", "rank")


household_MR <- function(household_GWAS_result, i, exposure_sex, grouping_var) {
  #household_GWAS_result <- "gwas_age_bins"
  #grouping_var <- "age_even_bins"

  trait_info <-  shiny_data[[i]][["trait_info"]]
  trait_description <- as.character(trait_info["description",1])
  phesant_ID <- as.character(trait_info["phes_ID",1])
  trait_ID <- as.character(trait_info["trait_ID",1])

  if(exposure_sex=="male"){outcome_sex="female"}
  if(exposure_sex=="female"){outcome_sex="male"}
  if(exposure_sex=="male"){index="HOUSEHOLD_MEMBER1"}
  if(exposure_sex=="female"){index="HOUSEHOLD_MEMBER2"}
  opp_index <- ifelse(index=="HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2", "HOUSEHOLD_MEMBER1")
  ####### load IV list
  IV_data <- shiny_data[[i]][[paste0(exposure_sex,"_IVs")]]


  colnames(IV_data) <- c("SNP", "chr", "beta", "se", "pval", "other_allele", "effect_allele","eaf","samplesize")
  data_IV_format <- format_data(IV_data, type="exposure")
  outcome <- "phenotype"
  j <- which(models_to_run[["exposure_sex"]]==exposure_sex & models_to_run[["trait_ID"]]==trait_ID) # j represents the model
  outcome_gwas <- get(household_GWAS_result)[[j]]$outcome_gwas_out
  household_intervals <- levels(outcome_gwas[[grouping_var]])

  outcome_gwas <- outcome_gwas %>% mutate_if(is.factor,as.character) %>%
    mutate_at(vars(geno_index_beta, geno_index_se, geno_index_pval, AF), as.numeric)
  outcome_dat_full <- format_data(outcome_gwas, type="outcome",
                             snp_col = "SNP",
                             beta_col = paste0("geno_index_beta"),
                             se_col = paste0("geno_index_se"),
                             effect_allele_col = "allele1",
                             other_allele_col = "allele0",
                             pval_col = paste0("geno_index_pval"),
                             eaf_col = "AF",
                             phenotype_col = grouping_var
  )
  harmonise_dat_full <- harmonise_data(
    exposure_dat = data_IV_format,
    outcome_dat = outcome_dat_full, action=1
  )

  full_MR_summary <-  numeric()
  bin_summary <- numeric()

  household_intervals_num <- str_first_number(household_intervals)

  household_intervals_num[which(is.na(household_intervals_num))] <- 1e3
  for(bin in household_intervals[order(household_intervals_num)]){


    gwas_sub <- outcome_gwas[which(outcome_gwas[[grouping_var]]==bin),] %>% mutate_if(is.factor,as.character) %>%
      mutate_at(vars(geno_index_beta, geno_index_se, geno_index_pval, AF), as.numeric)

    n_bin <- max(gwas_sub$n)
    outcome_dat <- format_data(gwas_sub, type="outcome",
                               snp_col = "SNP",
                               beta_col = paste0("geno_index_beta"),
                               se_col = paste0("geno_index_se"),
                               effect_allele_col = "allele1",
                               other_allele_col = "allele0",
                               pval_col = paste0("geno_index_pval"),
                               eaf_col = "AF"
    )

    harmonise_dat <- harmonise_data(
      exposure_dat = data_IV_format,
      outcome_dat = outcome_dat, action=1
    )


    original_MR <- mr(harmonise_dat, method=MR_method_list)
    if(dim(original_MR)[1]!=0){
      MR_res <- include_MR_NA(original_MR)
      MR_ivw_row <- which(MR_res[,"method"]=="Inverse variance weighted")
      MR_wald_row <- which(MR_res[,"method"]=="Wald ratio")
      #MR_dir <- paste0(pheno_dir, "/household_MR/exposure_", exposure_sex)
      #write.csv(MR_res, paste0(MR_dir, "/", phenotype_description, "_", exposure_sex,"-",outcome_sex, "_MR_bin_",bin, ".csv"), row.names=F)

      bin_result <- c(bin, outcome, n_bin, exposure_sex, outcome_sex, MR_res[MR_ivw_row,"b"],MR_res[MR_ivw_row,"se"], MR_res[MR_ivw_row,"pval"])
    } else bin_result <- c(bin, outcome, n_bin, exposure_sex, outcome_sex, NA,NA, NA)
    bin_summary <- rbind(bin_summary, bin_result)



    if(bin=="all"){

    het_test <- mr_heterogeneity(harmonise_dat, method_list=c("mr_egger_regression", "mr_ivw"))
    egger_test <- mr_pleiotropy_test(harmonise_dat)
    leave1out_test <- mr_leaveoneout(harmonise_dat)

    #write.csv(het_test, paste0(MR_dir, "/", phenotype_description, "_", exposure_sex,"-",outcome_sex, "_MR-het.csv"), row.names=F)
    #write.csv(egger_test, paste0(MR_dir, "/", phenotype_description, "_", exposure_sex,"-",outcome_sex, "_MR-Egger.csv"), row.names=F)
    #write.csv(leave1out_test, paste0(MR_dir, "/", phenotype_description, "_", exposure_sex,"-",outcome_sex, "_MR-leave_out.csv"), row.names=F)

    gwas_sub <- outcome_gwas[which(outcome_gwas[[grouping_var]]==bin),] %>% mutate_if(is.factor,as.character) %>%
      mutate_at(vars(geno_index_beta, geno_index_se, geno_index_pval, AF), as.numeric)
    ngwas <- max(gwas_sub$n, na.rm=T)
    nsnps <- dim(harmonise_dat)[1]
    n_neale <- max(data_IV_format$samplesize.exposure, na.rm=T)
    temp_summary <- summarize_mr_result (paste0(trait_ID,"_INDEX"), paste0(trait_ID,"_HOUSEHOLD"), exposure_sex, outcome_sex, nsnps,MR_res, het_test, egger_test,n_neale, ngwas)

    ## Test for reverse causality
    harmonise_dat_sensitivity <- harmonise_dat[which(harmonise_dat[["pval.exposure"]] < harmonise_dat[["pval.outcome"]]),]
    MR_res <- include_MR_NA(mr(harmonise_dat_sensitivity, method=MR_method_list))
    #write.csv(MR_res, paste0(MR_dir, "/", phenotype_description, "_", exposure_sex,"-",outcome_sex, "_MR-sensitivity.csv"), row.names=F)

    MR_ivw_row <- which(MR_res[,"method"]=="Inverse variance weighted")
    MR_wald_row <- which(MR_res[,"method"]=="Wald ratio")
    nsnps_sensitivity <- dim(harmonise_dat_sensitivity)[1]
    correct_row <- ifelse(nsnps_sensitivity==1, MR_wald_row, MR_ivw_row)
    sensitivity_result <- cbind(nsnps_sensitivity, make_beta_95ci(MR_res[correct_row,"b"],MR_res[correct_row,"se"]),pretty_round(MR_res[correct_row,"pval"]))
    temp_summary <- cbind(temp_summary, sensitivity_result)

    ## Make MR plot
    #output_figure_dir <- paste0(project_dir, "/output/figures/traitMR/",trait_ID)
    #pdf(file=paste0(output_figure_dir,"/", trait_ID, "_", exposure_sex, "-",outcome_sex, "_household_MR", ".pdf"))
    mr_title <- bquote(atop(.(paste0("Estimate of the assortative mating effect of ")),
                       italic(.(trait_description)) ~ .(paste0(' for ', exposure_sex, "s on ", outcome_sex,"s"))))
    mr_plot <- mr_scatter_plot_custom(original_MR, harmonise_dat, mr_title, exposure_sex, outcome_sex)
    #dev.off()

    MR_summary <- cbind(outcome, temp_summary)
    full_MR_summary <- rbind(full_MR_summary, MR_summary)
  }
  }

  num_cols <- length(colnames(full_MR_summary))

  colnames(full_MR_summary)[num_cols-1] <- "IVW/Wald_summary_sensitivity"
  colnames(full_MR_summary)[num_cols] <- "IVW/Wald_pval_sensitivity"
  colnames(full_MR_summary)[num_cols-2] <- "N_snps_sensitivity" #test for reverse causation
  #output_table_dir <- paste0(project_dir, "/output/tables/traitMR/", trait_ID)
  #write.csv(full_MR_summary, paste0(output_table_dir, "/", trait_ID, "_household_MR.csv"), row.names=F, quote=T)

  colnames(bin_summary) <- c("bin","outcome","n", "exposure_sex","outcome_sex", "IVW_beta", "IVW_se", "IVW_pval")
  #write.csv(bin_summary, paste0(pheno_dir, "/household_MR/", phenotype_description,"_", exposure_sex, "-",outcome_sex, "_household_MR_bin.csv"), row.names=F, quote=T)

  out <- list(harmonise_dat_full = harmonise_dat_full, bin_summary = bin_summary, full_MR_summary = full_MR_summary, mr_plot = mr_plot, leave1out_test = leave1out_test)

}

bin_plot <- function(household_mr_output, grouping_var, i){

  harmonise_dat_full <- household_mr_output$harmonise_dat_full %>% dplyr::filter(outcome != "all") %>% mutate_at('outcome', as.factor)

  outcome_levels <- levels(harmonise_dat_full$outcome)
  outcome_levels_num <- str_first_number(outcome_levels)

  harmonise_dat_full$outcome <- factor(harmonise_dat_full$outcome, levels = outcome_levels[order(outcome_levels_num)])

  bin_summary <- household_mr_output$bin_summary %>% as_tibble %>% dplyr::filter(bin != "all") %>% mutate_at('bin', as.factor)
  bin_summary$bin <- factor(bin_summary$bin, levels = outcome_levels[order(outcome_levels_num)])


  trait_info <-  shiny_data[[i]][["trait_info"]]
  trait_description <- as.character(trait_info["description",1])

  legend_title <- ifelse(grouping_var=="age_even_bins", "Median age \n of couples", "Time together \n in same household")
  bin_summary <- as.data.frame(bin_summary) %>% mutate_if(is.factor,as.character) %>%
    mutate_at(vars(IVW_beta, IVW_se, IVW_pval), as.numeric)

  plot <- ggplot2::ggplot(data = harmonise_dat_full, ggplot2::aes(x = beta.exposure,
                                         y = beta.outcome)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome),
                                                                                     colour = "grey", width = 0) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), colour = "grey", height = 0) +
    ggplot2::geom_point(ggplot2::aes(colour = factor(outcome))) +
    ggplot2::geom_abline(data = bin_summary, ggplot2::aes(intercept = 0,
                                                    slope = IVW_beta, colour = factor(bin)), show.legend = TRUE) +
    ggplot2::scale_colour_manual(values = brewer.pal(n = 9, name = "Blues")[-c(1:2)]) +
    theme_minimal() +
    ggplot2::labs(colour = legend_title, x = paste("SNP effect on", trait_description),
                  y = paste("SNP effect on", "partner")) +
    ggplot2::theme(legend.position = "right", legend.direction = "vertical") +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1))
  return(plot)
}

mr_scatter_plot_custom <-  function (mr_results, dat, mr_title, exposure_sex, outcome_sex)
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"),
                       function(d) {
                         d <- plyr::mutate(d)
                         if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
                           return(blank_plot("Insufficient number of SNPs"))
                         }
                         d <- subset(d, mr_keep)
                         index <- d$beta.exposure < 0
                         d$beta.exposure[index] <- d$beta.exposure[index] *
                           -1
                         d$beta.outcome[index] <- d$beta.outcome[index] *
                           -1
                         mrres <- subset(mr_results, id.exposure == d$id.exposure[1] &
                                           id.outcome == d$id.outcome[1])
                         mrres$a <- 0
                         if ("MR Egger" %in% mrres$method) {
                           temp <- mr_egger_regression(d$beta.exposure,
                                                       d$beta.outcome, d$se.exposure, d$se.outcome,
                                                       default_parameters())
                           mrres$a[mrres$method == "MR Egger"] <- temp$b_i
                         }
                         if ("MR Egger (bootstrap)" %in% mrres$method) {
                           temp <- mr_egger_regression_bootstrap(d$beta.exposure,
                                                                 d$beta.outcome, d$se.exposure, d$se.outcome,
                                                                 default_parameters())
                           mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
                         }
                         ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure,
                                                                y = beta.outcome)) +
                           ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome -se.outcome, ymax = beta.outcome + se.outcome),
                                                                                                            colour = "grey", width = 0) +
                           ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), colour = "grey", height = 0) +
                           ggplot2::geom_point() + #ggplot2::aes(text = paste("SNP:", SNP))) +
                           ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, slope = b, colour = method), show.legend = TRUE) +
                           theme_minimal() + labs(title=mr_title) + theme(plot.title = element_text(hjust = 0.5)) +
                           ggplot2::scale_colour_manual(values = c("#a6cee3",
                                                                   "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
                                                                   "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                                                                   "#6a3d9a", "#ffff99", "#b15928")) + ggplot2::labs(colour = "MR Test",
                                                                                                                     x = paste0("SNP effect on ", exposure_sex, "s (general population)"), y = paste("SNP effect on",
                                                                                                                                                                          outcome_sex, "partner")) + ggplot2::theme(legend.position = "bottom",
                                                                                                                                                                                                          legend.direction = "horizontal") + ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
                           theme(legend.title = element_blank())
                       })
  mrres
}


include_MR_NA <- function(MR_res){
  out <- MR_res
  used_methods <- MR_res[,"method"]
  method_col <- which(colnames(MR_res)=="method")
  total_col <- dim(MR_res)[2]
  attempted_methods <- levels(used_methods)[which(!levels(used_methods) %in%  used_methods)]
  for(i in attempted_methods)
  {
    empty_row <- c(rep(NA,method_col-1), i, rep(NA,total_col-method_col))
    out <- rbind(out, empty_row)
  }
  return(out)
}

summarize_mr_result <- function(exposure, outcome, exposure_sex, outcome_sex, nsnps, MR_res, het_test, egger_test,n_exposure,n_outcome){
  out <- cbind(exposure, outcome, exposure_sex, outcome_sex, nsnps, n_exposure,n_outcome)
  if(nsnps==1)
  {
    egger_pval <- NA
    het_IVW_pval <- NA
    for(method in c("Wald ratio","MR Egger"))
    {
      row <- which(MR_res[,"method"]==method)
      if(any(is.na(MR_res[row,c("b","se","pval")])))
      {
        out <- cbind(out, NA, NA)
      } else {
        b.ci <- make_beta_95ci(MR_res[row,"b"],MR_res[row,"se"])
        pval <- pretty_round(MR_res[row,"pval"])
        out <- cbind(out, b.ci, pval)
      }
    }
    out <- cbind(out, egger_pval,het_IVW_pval)
  }
  if(nsnps>1)
  {

    for(method in c("Inverse variance weighted","MR Egger"))
    {
      row <- which(MR_res[,"method"]==method)
      b.ci <- make_beta_95ci(MR_res[row,"b"],MR_res[row,"se"])
      pval <- pretty_round(MR_res[row,"pval"])
      out <- cbind(out, b.ci, pval)
    }
    het_IVW_pval <- pretty_round(het_test[which(het_test[,"method"]=="Inverse variance weighted"),"Q_pval"])
    egger_pval <- pretty_round(egger_test[1,"pval"])
    out <- cbind(out, egger_pval,het_IVW_pval)
  }
  colnames(out) <- c("exposure", "outcome", "exposure_sex","outcome_sex","N_snps", "N_exposure_GWAS", "N_outcome_GWAS",
                     "IVW/Wald_summary", "IVW/Wald_pval","Egger_summary", "Egger_pval", "Egger_int_pval", "Het_IVW_Qpval")
  return(out)
}


make_beta_95ci <- function(beta, se){
  beta <- as.numeric(beta)
  se <- as.numeric(se)
  beta_round <- pretty_round(beta)
  lower_ci <- pretty_round(beta-(1.96*se))
  upper_ci <- pretty_round(beta+(1.96*se))
  out <- paste0(beta_round," (",lower_ci,", ", upper_ci, ")")
  return(out)
}


pretty_round <- function(x){
  x <- as.numeric(x)
  if(!is.na(x))
  {
    if(abs(x)>1) {out<- sprintf("%.2f", round(x,2))}
    if(abs(x)<=1 & abs(x)>=0.001) {
      om = floor(log10(abs(x)))
      dp = 2-om-1
      out<- sprintf(paste("%.",dp,"f", sep=""), signif(x,2))}
    if(abs(x)<0.001) {out <- sprintf("%.1e", signif(x, 2))}
  } else out <- NA
  return(out)
}


