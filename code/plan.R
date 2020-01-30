
# Workflow Plans
#

## Description: creates pheno file and household pairs file for downstream analyses

tidy <- drake_plan(

  ##################################
  ### MAKE HOUSEHOLD PAIRS FILE ####
  ##################################
  household_info = read.csv(file_in(!!paste0(UKBB_dir,"/pheno/ukb6881.csv")), header=T),
  relatives = read.table(file_in(!!paste0(UKBB_dir,"/geno/","ukb1638_rel_s488366.dat")), header=T),
  phesant_directory = read.table(file_in(!!paste0(UKBB_processed,"/PHESANT/","PHESANT_file_directory.txt")), header=T),
  phesant_files = file_in(!!as.character(phesant_file_list)),
  hh_pairs = pairs_only(household_info),
  hh_pairs_kin = find_kinship(hh_pairs, relatives),

  household_relationships = fread(file_in(!!relations_file), header=T, select=c("userId","sex",!!relatedness_field), data.table=F),
  hh_pairs_filter = filter_pairs(hh_pairs_kin,household_relationships,!!relatedness_field),
  #create directory
  write_hh_pairs = write.csv(hh_pairs_filter, file_out("analysis/data_setup/household_pairs.csv"),row.names=F, quote=T ),
  ################################
  ######## CREATE PHENO FILE #####
  ################################
  pheno = fread(file_in(!!first_phesant_file), header=T, select=c("userId","age"), data.table=F),
  sqc = fread(file_in(!!paste0(UKBB_dir,"/geno/ukb_sqc_v2.txt")), header=F, data.table=F),
  fam = fread(file_in(!!paste0(UKBB_dir,"/plink/_001_ukb_cal_chr9_v2.fam")), header=F,data.table=F),
  sqc_munge = munge_sqc(sqc,fam),
  model_adjustments = add_pcs(hh_pairs_filter,pheno,sqc_munge),

  write_houshold1 = write.csv(model_adjustments[[1]], file_out("analysis/data_setup/HOUSEHOLD_MEMBER1_pheno_model_adjustments.csv"),row.names=F, quote=T ),
  write_houshold2 = write.csv(model_adjustments[[2]], file_out("analysis/data_setup/HOUSEHOLD_MEMBER2_pheno_model_adjustments.csv"),row.names=F, quote=T ),
  joint_model_adjustments = merge(model_adjustments[[1]], model_adjustments[[2]], by=c("HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2", "HOUSE_ID","kinship","HOUSEHOLD_MEMBER1_sex","HOUSEHOLD_MEMBER2_sex")),
  write_joint = write.csv(joint_model_adjustments, file_out("analysis/data_setup/JOINT_pheno_model_adjustments.csv"), row.names=F, quote=T ),
  #####################################
  #### ADD TIME AT HOUSEHOLD ##########
  #####################################
  time_at_household = fread(file_in(!!time_at_address_file),  select=c("userId",time_at_address_field), data.table=F),
  time_at_household_raw = fread(file_in(!!time_at_address_raw_file),  data.table=F),
  household_time = calc_time_together(joint_model_adjustments,time_at_household,time_at_household_raw),
  household_intervals = {
    start <- 0
    end <- 50
    intervals <- c(seq(!!time_together_min,!!time_together_max, by=!!time_together_interval))
    labels <- interval_labels(intervals)
  },
  household_time_munge = {

    munge_household_time(household_time, joint_model_adjustments, !!time_together_min,!!time_together_max, !!time_together_interval,
    !!age_min, !!age_max, !!age_interval, !!num_household_bins)
  },
  write_household_time = write.csv(household_time_munge, file_out("analysis/data_setup/household_time.csv"), row.names=F, quote=T )

)

filter_traits <- drake_plan(

  ## old A1
  UKBB_directory = read.csv(file_in(!!paste0(UKBB_processed,"/UKBB_pheno_directory.csv")), header=T),
  trait_corrs = {
    phesant_files
    compute_trait_corr(phesant_directory,UKBB_directory,hh_pairs_filter)
  },
  write_traits_corr = write.csv(trait_corrs, file_out( "output/tables/1.household_correlations.csv"), row.names=F),

  ## old A2
  Neale_SGG_dir =  read.csv(!!paste0(Neale_summary_dir,"/Neale_SGG_directory.csv"), header=T), ## not tracking this as input because it changes overtime
  Neale_SGG_dir_filt = SGG_link_with_Neale(Neale_SGG_dir),
  traits_corr2 = filter_by_corr(trait_corrs,Neale_SGG_dir_filt,!!household_correlation_threshold),
  Neale_to_process = organize_Neale(traits_corr2),
  write_define_cats = write.csv(Neale_to_process$define_cats, file_out("output/tables/define_Neale_categories.csv"), row.names=F),
  write_download_list = write.csv(Neale_to_process$download_rest, file_out("analysis/download_Neale_list.csv"), row.names=F),

  ## old A3
  filled_cats = read.csv(file_in("output/tables/define_Neale_categories_filled.csv"), check.names=F),
  Neale_manifest = read_tsv(file_in(!!paste0(Neale_output_path,"/",Neale_manifest_file_name)), col_types = cols()),
  traits_corr2_filled = download_Neale(filled_cats,Neale_to_process$download_rest,traits_corr2,
    file_in(!!paste0(Neale_output_path,"/",Neale_manifest_file_name))),
  #write_traits_corr2 = write.csv(traits_corr2_filled, "output/tables/2.household_correlations.corr_filter.csv", row.names=F),

  ## old A4
  run_process_Neale = {
    traits_corr2_filled
    processx::run(command = "sbatch", c(file_in("code/process_Neale.sh")))},
  # this script wil update the `Neale_SGG_directory` which will make this file outdated

  ## old A5
  traits_corr2_update = {
    run_process_Neale
    stats1 <- 2
    stats2 <- 2
    while ( (stats1 > 1 | stats2 > 1)){
      stats1 <- suppressWarnings(system(paste("squeue -n", "process_Neale"), intern = TRUE))
      stats2 <- suppressWarnings(system(paste("squeue -n", "clump_Neale_IVs"), intern = TRUE))
      print("Still running...")
      Sys.sleep(60)
    }

     update_download_info(traits_corr2_filled,Neale_SGG_dir)
  },
  traits_to_count_IVs = tibble(
      i = 1:dim(traits_corr2_update[which(traits_corr2_update[["Neale_file_sex"]]=="both"),])[1],
    ),
  clump_dir = target(!!paste0(Neale_summary_dir,"/IVs/clump/" )),
    #trigger = trigger(change = file.mtime(!!paste0(Neale_summary_dir,"/IVs/clump/" )))),
  Neale_files_dir = target(!!Neale_output_path),
    #trigger = trigger(change = file.mtime(!!Neale_output_path))),
  IV_list = target(
    {
      clump_dir
      #files <- file_in_out(traits_corr2_update,traits_to_count_IVs$i, Neale_manifest, !!IV_threshold)
      get_IV_list(traits_corr2_update,traits_to_count_IVs$i,Neale_manifest,!!IV_threshold)
    }, dynamic = map(traits_to_count_IVs)
  ),

  IV_list_out = write_IV_list(traits_corr2_update, traits_to_count_IVs, IV_list,
    !!IV_threshold, file_out("analysis/data_setup/IV_lists/")),
  ## old A6
  IV_summary = target(
    {
      IV_list_out
      Neale_files_dir
      summarize_IVs(traits_corr2_update,traits_to_count_IVs$i,Neale_manifest)
    }, dynamic = map(traits_to_count_IVs)
  ), ## the output of IV_summary is a matrix
  traits_corr3 = IV_filter(traits_corr2_update, IV_summary, !!num_IVs_threshold),
  write_traits_corr2 = write.csv(traits_corr3$non_filtered,file_out("output/tables/2.household_correlations.corr_filter.csv"), row.names=F),

  ## old A7
  variant_data_reduced = target(reduce_variant_data(traits_corr3$to_run,file_in(!!variant_file_full_name)), hpc = FALSE),
  traits_to_calc_het = tibble(
      i = 1:dim(traits_corr3$to_run)[1],
    ),
  sex_het_summary = target({
    trait_ID <- as.character(traits_corr3$to_run[traits_to_calc_het$i,"Neale_pheno_ID"])
    calc_sex_het(traits_corr3$to_run,traits_to_calc_het$i,variant_data_reduced,Neale_manifest)
  }, dynamic = map(traits_to_calc_het)
  ),

  IV_info_out = write_IV_info(sex_het_summary, traits_to_calc_het, traits_corr3,
    file_out("analysis/data_setup/IV_info/")),

  sex_het_out = write_sex_het(sex_het_summary, traits_to_calc_het, traits_corr3, file_out("analysis/data_setup/sex_heterogeneity/")),
  traits_corr4 = sex_het_filter(traits_corr3$to_run, sex_het_summary, traits_to_calc_het, !!num_IVs_threshold),
  write_traits_corr3 = write.csv(traits_corr4$non_filtered, file_out("output/tables/3.household_correlations.numIVs_filter.csv"), row.names=F),
  write_traits_corr4 = write.csv(traits_corr4$to_run, file_out("output/tables/4.household_correlations.sexhet_filter.csv"), row.names=F),

  ## old A8
  valid_GRS_summary = check_valid_GRS_input(traits_corr4$to_run,Neale_manifest),
  traits_corr5 = valid_GRS_filter(traits_corr4$to_run,valid_GRS_summary),
  write_traits_corr5 = write.csv(traits_corr5, file_out("output/tables/5.household_correlations.baseGRS_filter.csv"), row.names=F),

  # filter for continuous
  traits_corr6 = traits_corr5 %>% filter(variable_type=="ordinal" | variable_type=="continuous_irnt"),
  write_traits_corr6 = write.csv(traits_corr6, file_out("output/tables/6.household_correlations.nonbinary_filter.csv"), row.names=F)

)

pre_pipeline <- bind_plans(tidy,filter_traits)

# config <- drake_config(pre_pipeline)
# vis_drake_graph(config)


pipeline <- drake_plan(

  create_final_filter = {
    file_in("output/tables/6.household_correlations.nonbinary_filter.csv")
    system("cp output/tables/6.household_correlations.nonbinary_filter.csv output/tables/household_correlations.final_filter.csv")
    file_out("output/tables/household_correlations.final_filter.csv")},
  traits = {
    create_final_filter
    read.csv(file_in("output/tables/household_correlations.final_filter.csv"), header=T)},
  traits_to_run = tibble(
      i = 1:dim(traits)[1],
    ),
  sample_file = fread(file_in(!!paste0(UKBB_dir, "/imp/ukb1638_imp_chr1_v2_s487398.sample")), skip=2, header=F,data.table=F),

  ## data prep
  data_prep = target(
    {
      phesant_files
      prep_data(traits,traits_to_run$i,phesant_directory,!!GRS_thresholds,Neale_manifest,sqc,fam)
    }, dynamic = map(traits_to_run)

  ),

  data_prep_out = target({
    data_prep
    write_data_prep(traits, traits_to_run, file_out("analysis/traitMR/trait_info"), file_out("analysis/traitMR/pheno_files/phesant"))
  }, hpc = FALSE),


  ## old create summary stats
  summ_stats_create = target(
    {
      data_prep_out
      trait_ID <- as.character(traits[traits_to_run$i,"Neale_pheno_ID"]) ## this is the Neale_id, used to be pheno_description
      file_in("analysis/traitMR/trait_info", "analysis/data_setup/sex_heterogeneity", "analysis/data_setup/IV_lists")
      create_summary_stats(traits,traits_to_run$i,phesant_directory,!!GRS_thresholds,Neale_manifest,variant_data_reduced)
    }, dynamic = map(traits_to_run)
  ),

  summ_stats_out = target({
    summ_stats_create
    write_summ_stats(summ_stats_create, traits, traits_to_run, !!GRS_thresholds, file_out("analysis/traitMR/IVs"),file_out("analysis/traitMR/GRS"))
  },hpc = FALSE),

  imp_dir = target(!!paste0(UKBB_dir, "/imp"),
    trigger = trigger(change = file.mtime(!!!paste0(UKBB_dir, "/imp")))),

  models_to_run = {
    file_in("analysis/traitMR/trait_info")
    summ_stats_out
    define_models(traits)
  },

  gwas_tt_bins = target({
    summ_stats_create
    file_in("analysis/traitMR/trait_info")
    file_in("analysis/traitMR/IVs")
    imp_dir
    household_GWAS(models_to_run$i,models_to_run$trait_ID, models_to_run$exposure_sex,
      (models_to_run$phenotype_file),models_to_run$phenotype_col,models_to_run$phenotype_description,
      (models_to_run$IV_file),sample_file,pheno_cov=joint_model_adjustments,
      grouping_var="time_together_even_bins",household_time_munge,(models_to_run$gwas_outcome_file))

    }, dynamic = map(models_to_run)
  ),

  gwas_age_bins = target({
    summ_stats_create
    file_in("analysis/traitMR/trait_info")
    file_in("analysis/traitMR/IVs")
    imp_dir
    household_GWAS(models_to_run$i,models_to_run$trait_ID, models_to_run$exposure_sex,
      (models_to_run$phenotype_file),models_to_run$phenotype_col,models_to_run$phenotype_description,
      (models_to_run$IV_file),sample_file,pheno_cov=joint_model_adjustments,
      grouping_var="age_even_bins",household_time_munge,(models_to_run$gwas_outcome_file))

    }, dynamic = map(models_to_run)
  ),

  mr_age_bins = target({

    file_in("analysis/traitMR/trait_info")
    file_in("analysis/traitMR/IVs")

    household_MR(household_GWAS_result = gwas_age_bins, models_to_run$i, models_to_run$trait_ID,
      models_to_run$exposure_sex,grouping_var="age_even_bins", models_to_run$IV_file)
    }, dynamic = map(gwas_age_bins, models_to_run)
  ),

  mr_tt_bins = target({
    file_in("analysis/traitMR/trait_info")
    file_in("analysis/traitMR/IVs")
    household_MR(household_GWAS_result = gwas_tt_bins, models_to_run$i, models_to_run$trait_ID,
      models_to_run$exposure_sex, grouping_var="time_together_even_bins", models_to_run$IV_file)
    }, dynamic = map(gwas_tt_bins, models_to_run)
  ),

  shiny_data = target({
    file_in("analysis/traitMR/trait_info")
    file_in("analysis/traitMR/IVs")
    read_shiny_data(traits,traits_to_run$i)

  }, dynamic = map(traits_to_run)),

  shiny_data_out = target({
    loadd(traits)
    loadd(shiny_data)
    loadd(models_to_run)
    loadd(mr_age_bins)
    loadd(mr_tt_bins)
    save(traits, shiny_data, models_to_run, gwas_age_bins, gwas_tt_bins, file = file_out("code/shiny/data.RData"))
  }, hpc = FALSE)


)



full_proxymr <- bind_plans(pre_pipeline,pipeline)


run_mr <- drake_plan(

  gwas_age_bins = readd(gwas_age_bins),
  gwas_tt_bins = readd(gwas_tt_bins),
  models_to_run = readd(models_to_run),

  mr_age_bins = target({

    file_in("analysis/traitMR/trait_info")
    file_in("analysis/traitMR/IVs")

    household_MR(household_GWAS_result = gwas_age_bins, models_to_run$i, models_to_run$trait_ID,
      models_to_run$exposure_sex,grouping_var="age_even_bins", models_to_run$IV_file)
    }, dynamic = map(gwas_age_bins, models_to_run)
  ),

  mr_tt_bins = target({
    file_in("analysis/traitMR/trait_info")
    file_in("analysis/traitMR/IVs")
    household_MR(household_GWAS_result = gwas_tt_bins, models_to_run$i, models_to_run$trait_ID,
      models_to_run$exposure_sex, grouping_var="time_together_even_bins", models_to_run$IV_file)
    }, dynamic = map(gwas_tt_bins, models_to_run)
  ),


)
