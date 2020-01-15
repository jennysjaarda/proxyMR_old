
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
    #unique_phes_ids <- unique(phesant_directory[,2])
  UKBB_directory = read.csv(file_in(!!paste0(UKBB_processed,"/UKBB_pheno_directory.csv")), header=T),
  #pairs_filter = hh_pairs_filter
  trait_corrs = {

    file_in(unique(phesant_directory$File))
    compute_trait_corr(phesant_directory,UKBB_directory,hh_pairs_filter)
  },
  write_traits_corr = write.csv(trait_corrs, file_out( "output/tables/1.household_correlations.csv"), row.names=F),

  ## old A2 - to modify
  Neale_SGG_dir =  read.csv(file_in(!!paste0(Neale_summary_dir,"/Neale_SGG_directory.csv")), header=T),
  Neale_SGG_dir_filt = SGG_link_with_Neale(Neale_SGG_dir),
  traits_corr2 = filter_by_corr(trait_corrs,Neale_SGG_dir_filt,!!household_correlation_threshold),
  Neale_to_process = organize_Neale(traits_corr2),
  write_define_cats = write.csv(Neale_to_process$define_cats, file_out("output/tables/define_Neale_categories.csv"), row.names=F),
  write_download_list = write.csv(Neale_to_process$download_rest, file_out("analysis/download_Neale_list.csv"), row.names=F),

  ## old A3 - to modify
  filled_cats = read.csv(file_in("output/tables/define_Neale_categories_filled.csv"), check.names=F),
  traits_corr2_filled = download_Neale(filled_cats,Neale_to_process$download_rest,traits_corr2),
  write_traits_corr2 = write.csv(traits_corr2_filled, "output/tables/2.household_correlations.corr_filter.csv", row.names=F),
  run_process_Neale = processx::run(command = "sbatch", c("code/process_Neale.sh")),
  IV_list = { while (lenght(.....))
    


  }
  ## need to find a way to track this job
  stat <- suppressWarnings(system(paste("squeue -n", "clump_Neale_IVs")))
#check if length(stat) >1

stat <- suppressWarnings(system(paste("squeue -n", slr_job$jobname),
        intern = TRUE))
    if (length(stat) > 1) {
        cat(paste(c("Job running or in queue. Status:", stat),
            collapse = "\n"))


)





  ## old A4 - to modify: add A5 sbatch script
  processx::run(sbatch "scripts/A4.process_Neale.sh"),

  ## old A6

)

join_ <- bind_plans(tidy,filter_traits)
make(join_)

#vis_drake_graph(drake_config(filter_traits))
