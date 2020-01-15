
# Workflow Plans
#

## Description: creates pheno file and household pairs file for downstream analyses

tidy <- drake_plan(

  ##################################
  ### MAKE HOUSEHOLD PAIRS FILE ####
  ##################################
  household_info = read.csv(file_in(!!paste0(UKBB_dir,"/pheno/ukb6881.csv")), header=T),
  relatives = read.table(!!paste0(UKBB_dir,"/geno/","ukb1638_rel_s488366.dat"), header=T),
  phesant_directory = read.table(file_in(!!paste0(UKBB_processed,"/PHESANT/","PHESANT_file_directory.txt")), header=T),
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
  joint_model_adjustments = merge(model_adjustments[[1]], model_adjustments[[1]], by=c("HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2", "HOUSE_ID","kinship","HOUSEHOLD_MEMBER1_sex","HOUSEHOLD_MEMBER2_sex")),
  write_joint = write.csv(joint_model_adjustments, file_out("analysis/data_setup/JOINT_pheno_model_adjustments.csv"), row.names=F, quote=T ),
  #####################################
  #### ADD TIME AT HOUSEHOLD ##########
  #####################################
  time_at_household = fread(file_in(!!time_at_address_file),  select=c("userId",time_at_address_field), data.table=F),
  time_at_household_raw = fread(file_in(!!time_at_address_raw_file),  data.table=F),
  household_time = calc_time_together(joint_model_adjustments,time_at_household,time_at_household_raw),
  write_household_time = write.csv(household_time, file_out("analysis/data_setup/household_time.csv"), row.names=F, quote=T )

)

filter_traits <- drake_plan(

  ## old A1

  UKBB_directory = read.csv(file_in(!!paste0(UKBB_processed,"/UKBB_pheno_directory.csv")), header=T),

  trait_corrs = {

    file_in(unique(phesant_directory$File))
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
  reference_file_name = file_in(!!file.path(Neale_output_path,Neale_manifest_file_name)),
  reference_file = read_tsv(reference_file_name, col_types = cols()),
  traits_corr2_filled = download_Neale(filled_cats,Neale_to_process$download_rest,traits_corr2, reference_file_name),
  #write_traits_corr2 = write.csv(traits_corr2_filled, "output/tables/2.household_correlations.corr_filter.csv", row.names=F),

  ## old A4
  run_process_Neale = processx::run(command = "sbatch", c(file_in("code/process_Neale.sh"))),
  # this script wil update the `Neale_SGG_directory` which will make this file outdate

  ## old A5
  traits_corr2_update = {
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
  IV_list = target({

      #file_in(!!Neale_output_path)
      #file_in(!!paste0(Neale_summary_dir,"/IVs/clump/" ))
      get_IV_list(traits_corr2_update,traits_to_process$i,reference_file)},
      dynamic = map(traits_to_count_IVs)
    ),

  ## old A6
  IV_summary = target({

    #file_in(!!Neale_summary_dir)
    summarize_IVs(traits_corr2_update,traits_to_process$i,reference_file)},
    dynamic = map(traits_to_count_IVs)
  ), ## the output of IV_summary is a matrix

  traits_corr3 = IV_filter(traits_corr2_update, IV_summary, !!num_IVs_threshold),
  write_traits_corr2 = write.csv(traits_corr3$non_filtered,file_out("output/tables/2.household_correlations.corr_filter.csv"), row.names=F),

  ## old A7
  variant_data = fread(file_in(!!variant_file_full_name),data.table=F),
  variant_data_reduced = reduce_variant_data(variant_data)
  traits_to_calc_het = tibble(
      i = 1:dim(traits_corr3$to_run)[1],
    ),
  sex_het_summary = target({
    calc_sex_het(traits_corr3$to_run,i,variant_data_reduced,output_folder = file_out("analysis/data_setup/sex_heterogeneity/"))

  }, dynamic = map(traits_to_calc_het)
  ),
  traits_corr4 = sex_het_filter(traits_corr3$to_run,sex_het_summary,!!num_IVs_threshold),
  write_traits_corr3 = write.csv(traits_corr4$non_filtered, file_out("3.household_correlations.numIVs_filter.csv"), row.names=F),
  write_traits_corr4 = write.csv(traits_corr4$to_run, file_out("4.household_correlations.sexhet_filter.csv"), row.names=F),

  ## old A8
  valid_GRS_summary = check_valid_GRS_input(traits_corr4$to_run,reference_file)
  traits_corr5 = valid_GRS_filter(traits_corr4$to_run,valid_GRS_summary)
  write.csv(traits_corr5, file_out("output/tables/5.household_correlations.baseGRS_filter.csv"), row.names=F)

)

errorFile <- paste0("get_IV_list", "_errors.txt")
if(file.exists(errorFile)) {file.remove(errorFile)}

join_ <- bind_plans(tidy,filter_traits)
#make(join_)

make(join_,memory_strategy = "lookahead", garbage_collection = TRUE, parallelism = "clustermq",console_log_file = "pipeline_prep.out", jobs = 20, template = list(cpus = 1, partition = "sgg"))
make(join_,memory_strategy = "lookahead", garbage_collection = TRUE, parallelism = "future",console_log_file = "pipeline_prep.out", jobs = 20)

pipeline <- drake_plan(


)

#vis_drake_graph(drake_config(filter_traits))
