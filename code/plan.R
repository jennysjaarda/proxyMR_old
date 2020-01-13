
# Workflow Plans
#

filter_traits <- drake_plan(

  ## old A1
  phesant_directory = read.table(file_in(paste0(!!UKBB_processed,"/PHESANT/","PHESANT_file_directory.txt")), header=T)
  #unique_phes_ids <- unique(phesant_directory[,2])
  UKBB_directory = read.csv(paste0(file_in(!!UKBB_processed,"/UKBB_pheno_directory.csv")), header=T)
  pairs_filter = read.csv(trait_corr, file_in("analysis/data_setup/household_pairs.csv"),header=T)
  traits = compute_trait_corr(phesant_directory)
  trait_corr_out = write.csv(file_out( "output/tables/1.household_correlations.csv"), row.names=F)

  ## old A2 - to modify
  Neale_SGG_dir =  read.csv(file_in(paste0(!!Neale_summary_dir,"/Neale_SGG_directory.csv")), header=T)
  Neale_SGG_dir_filt2 = link_with_Neale(Neale_SGG_dir)
  traits_corr_filter = filter_by_corr(Neale_SGG_dir_filt2,traits)
  #write.csv(traits_corr_filter, "output/tables/2.household_correlations.corr_filter.csv", row.names=F)
  Neale_to_process = organize_Neale(traits_corr_filter)
  write.csv(Neale_to_process$define_cats2, paste0("output/tables/define_Neale_categories.csv"), row.names=F)
  write.csv(Neale_to_process$download_rest, paste0("pipeline/download_Neale_list.csv"), row.names=F)

  ## old A3 - to modify
  filled_cats <- read.csv("output/tables/define_Neale_categories_filled.csv", check.names=F)
  traits_corr_filter2 = download_Neale(filled_cats,Neale_to_process$download_rest,traits_corr_filter)
  write.csv(traits_corr_filter2, "output/tables/2.household_correlations.corr_filter.csv", row.names=F)

  ## old A4 - to modify: add A5 sbatch script
  processx::run(sbatch "scripts/A4.process_Neale.sh")

  ## old A6

)
