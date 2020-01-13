### File A1 ----

## Description: tests phenotypic correlations between household pairs that have been
## filtered for kinship < 0.05, and to be opposite-sex pairs
## run within project folder: /data/sgg2/jenny/projects/MR_Shared_Environment

compute_trait_corr <- function(phesant_directory){
  unique_phes_ids <- unique(phesant_directory[,2])
  trait_corr <- numeric()
  prev_file <- ""
  for (id in unique_phes_ids)
  {
      id <- as.character(id)
      phenotype_sub <- unlist(strsplit (id,"_"))[1]

      UKBB_directory[,"field.showcase"] <- as.character(UKBB_directory[,"field.showcase"])
      file_with_pheno <- which(UKBB_directory[,"field.showcase"]==id)
      file_name <- names(which.max(table(as.character(UKBB_directory[file_with_pheno,"file"]))))
      if(length(file_with_pheno)==0)
      {
        file_with_pheno <- which(UKBB_directory[,"field.showcase"]==phenotype_sub)
        file_name <- names(which.max(table(as.character(UKBB_directory[file_with_pheno,"file"]))))
      }
      description <- unlist(strsplit (as.character(UKBB_directory[file_with_pheno[1],"col.name"]),"_datacoding"))[1]
      current_file <- as.character(phesant_directory[which(phesant_directory[,2]==id),"File"])[1]

      #sgg_file <- sgg_paths[grepl(file_name,current_file)][1]
      #flow <- gsub(basename(sgg_file),"variable-flow-counts-all.txt", sgg_file)
      #file_dir <- gsub("/([^/]*)$" , "", sgg_file)
      #bin_num <- substrRight(file_dir,1)
      #log <- read.table(paste0(file_dir,"/out_bin",bin_num, "..log"))

      if(prev_file!=current_file)
      {tsv_data <- fread(current_file, header=TRUE, sep='\t',data.table=F)}
      r2 <- NA
      pairs_filter_copy <- pairs_filter
      pairs_filter_copy$trait1 <- tsv_data[[id]][match(pairs_filter_copy[["HOUSEHOLD_MEMBER1"]], tsv_data$userId)]
      pairs_filter_copy$trait2 <- tsv_data[[id]][match(pairs_filter_copy[["HOUSEHOLD_MEMBER2"]], tsv_data$userId)]
      complete_pairs <- pairs_filter_copy[which(complete.cases(pairs_filter_copy$trait1 ,pairs_filter_copy$trait2 )),]
      n_completed_pairs <- length(pairs_filter_copy[which(complete.cases(pairs_filter_copy$trait2, pairs_filter_copy$trait1)),"trait2"])
      if(n_completed_pairs>1)
      {
        sd1 <- sd(complete_pairs$trait1,na.rm=TRUE)
        sd2 <- sd(complete_pairs$trait2,na.rm=TRUE)
        if(sd1!=0 & sd2!=0)
          {
            r <- cor(pairs_filter_copy$trait1, pairs_filter_copy$trait2, method = c("pearson"),use = "pairwise.complete.obs")
            r2 <- r^2
          }
      }

      prev_file <- current_file
      trait_row <- cbind(id, phenotype_sub,description, n_completed_pairs,r2)
      trait_corr <- rbind(trait_corr, trait_row)
  }
  colnames(trait_corr) <- c("ID", "ID_sub", "description", "N_pairs","r2")
  cat(paste0("Household phenotypic correlations successfully computed, and saved to:\n",
  "'output/tables/1.household_correlations.csv'.\n"))
  return(trait_corr)
}


## File A2 ----

## Description: filter phenotype correlations for > specified threshold (in settings).
## Check which phenotypes have been downloaed, and create list of need to be downloaded
## and defined here: "output/tables/define_Neale_categories.csv".
## Manually categorize this list and save as: "output/tables/define_Neale_categories_filled.csv"
## run within project folder: /data/sgg2/jenny/projects/MR_Shared_Environment

link_with_Neale <- function(Neale_SGG_dir){

  Neale_SGG_dir_filt <- subset(Neale_SGG_dir, Neale_SGG_dir$phesant_processed=="YES")
  unique_Neale_ids <- unique((Neale_SGG_dir_filt[,"phenotype"]))

  both_sex_avail <- numeric()
  for(id in unique_Neale_ids)
  {
    id <- as.character(id)
    sex <- c(as.character(Neale_SGG_dir_filt[which(Neale_SGG_dir_filt[,"phenotype"]==id),"phenotype_sex"]))
    if("male" %in% sex & "female" %in% sex & "both" %in% sex)
    {
      both_sex_avail <- c(both_sex_avail,id)
    }
  }
  cat(paste0("There are: ", length(both_sex_avail), " phenotypes with Neale summary stats with joint
  and sex-specific data that have also been processed in the SGG database.\n\n")) ###1243 traits
  Neale_SGG_dir_filt2 <- subset(Neale_SGG_dir_filt, Neale_SGG_dir_filt$phenotype %in% both_sex_avail)

}

filter_by_corr <- function(trait_corr,Neale_SGG_dir_filt2){

  merge_temp <- merge(trait_corr, Neale_SGG_dir_filt2, by.x="ID", by.y="sgg_phesant_name", fill=T)
  merge_temp$r2 <- as.numeric(as.character(merge_temp$r2))
  corr_traits <- merge_temp[which(sqrt(merge_temp$r2) > household_correlation_threshold),]
  corr_traits <- corr_traits[,-which(colnames(corr_traits) %in% c("description.x"))]
  colnames(corr_traits) <- c("SGG_PHESANT_ID", "SGG_PHESANT_ID_sub","N_pairs","r2",
              "Neale_pheno_ID","Neale_pheno_ID_sub", "Neale_file_sex","description",
              "variable_type", "SGG_request", "PHESANT_processed", "SGG_location",
              "Neale_downloaded","category","define_category:T/F",
              "v2_exists", "v2_downloaded")
}

organize_Neale <- function(traits_corr_filter){
  corr_traits_sub <- subset(traits_corr_filter, traits_corr_filter$Neale_file_sex=="both")
  define_cats <- corr_traits_sub[,c("Neale_pheno_ID","Neale_pheno_ID_sub","description","variable_type","Neale_downloaded","category","define_category:T/F")]
  download_rest <- subset(define_cats,define_cats[["Neale_downloaded"]]=="PARTIALLY" & define_cats[["define_category:T/F"]]==FALSE)
  #write.csv(corr_traits, "pull_Neale_phenotypes.csv", row.names=F)
  date <- Sys.Date()
  define_cats2 <- subset(define_cats,define_cats[["define_category:T/F"]]==TRUE)

  if(dim(define_cats2)[1]==0 & dim(download_rest)[1]==0)
  {

      cat("All necessary Neale files have been downloaded, proceed in pipeline.\n\n")
  }
  if(dim(define_cats2)[1]!=0)
  {

      cat(paste0("Some Neale files need to be categorized, please categorize appropriately within file:
      'output/tables/define_Neale_categories.csv'.\n\n"))
  }

  if(dim(download_rest)[1]!=0)
  {

      cat(paste0("Some Neale files need to be downloaded, but have been categorized previously,
      likely they were updated based on Neale modifications to the PHESANT pipeline.
      See here: https://www.dropbox.com/s/ayto55r0z1eigwa/phenotype_update_note.txt?dl=0, for details.\n\n"))
  }

  cat(paste0("Finished linking phenotypic correlations file with Neale files.\n\n"))

  cat(paste0("Household phenotypic correlations have been filtered to include only phenotypes
  which have valid Neale summary stats in both sexes seperately and together,
  and have household correlations greater than specified r2 threshold (see 'scripts/settings.r' for chosen threshold),\n",
  "filtered phenotypic correlations have been saved to:
  'output/tables/2.household_correlations.corr_filter.csv'.\n\n"))

  return(list(download, define_cats))
}

### File A3 ----

## Description: download missing neale files using categories defined here:
## "output/tables/define_Neale_categories_filled.csv".
## Update categories in phenotypic correlations file: output/tables/pheno_household_filter_corr.csv
## run within project folder: /data/sgg2/jenny/projects/MR_Shared_Environment

download_Neale <- function(filled_cats,download_rest,traits_corr_filter){
  full_dl_list <- rbind(filled_cats, download_rest)
  select <- dplyr::select

  for(category
     in c("body", "brain", "diet", "disease", "disease_proxy", "lifestyle", "parental_pheno"))
  {
  IDs <- as.character(full_dl_list[which(full_dl_list[,"category"]==category),"Neale_pheno_ID"])
  download_neale_files( IDs,
                      category = category , reference_file=paste0(Neale_output_path, "/", Neale_manifest_file_name))
  }

  traits_corr_filter <-corr_traits
  levels(corr_traits$category) <- union (levels(corr_traits$category), levels(filled_cats$category))
  index <- which(corr_traits[["define_category:T/F"]]==TRUE)
  pheno_missing_cat <- corr_traits[index,"Neale_pheno_ID"]
  corr_traits[index, "category"] <- as.factor(filled_cats[["category"]][match(pheno_missing_cat, filled_cats$Neale_pheno_ID)])

  cat(paste0("Missing Neale files have been successfully downloaded and categorized,\n",
  "now they need to be processed (clumped for LD and extracted for p<0.1).\n"))


}
