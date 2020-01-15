pairs_only <- function(household_info){

  household_info <- household_info %>% mutate(group_count = ave(household_info$group, household_info$group,  FUN = length))
  #table(household_info$group_count) ## number of individuals in each group household sizes
  household_info <- household_info[which(household_info$group_count==2),]
  return(household_info)
}

find_kinship <- function(household_pairs, relatives){


  household_list <- unique(household_pairs$group)
  pairs <- numeric()
  for(house in household_list)
  {

    house_ids <- household_pairs[which(household_pairs$group==house),]
    for(i in 1:(dim(house_ids)[1]-1))
    {
      for(j in i:dim(house_ids)[1])
      {
        if (i==j) next
        else     (pairs <- rbind(pairs,cbind(house_ids[i,1],house_ids[j,1], house)))
      }
    }
  }

  colnames(pairs)[1] <- "HOUSEHOLD_MEMBER1"  ##Index
  colnames(pairs)[2] <- "HOUSEHOLD_MEMBER2"  ##Household member
  colnames(pairs)[3] <- "HOUSE_ID"  ##Household member

  pairs <- as.data.frame(pairs)
  pairs$kinship <- 0
  for (i in 1:dim(relatives)[1])
  {
    ind1 <- relatives[i,1]
    ind2 <- relatives[i,2]
    if(any((pairs[,1]==ind1 & pairs[,2]==ind2) | (pairs[,2]==ind1 & pairs[,1]==ind2)))
      {
        index <- which((pairs[,1]==ind1 & pairs[,2]==ind2) | (pairs[,2]==ind1 & pairs[,1]==ind2))
        pairs$kinship[index] <- relatives[i,5]
        print(i)
      }
  }
  return(pairs)
}




filter_pairs <- function(pairs,household_relationships,field){
  pairs$sex1 <- household_relationships[["sex"]][match(pairs[["HOUSEHOLD_MEMBER1"]], household_relationships$userId)]
  pairs$sex2 <- household_relationships[["sex"]][match(pairs[["HOUSEHOLD_MEMBER2"]], household_relationships$userId)]

  #merge(pairs, household_relationships, by.x=, by.y=)
  pairs_full_sex <- pairs[which(!is.na(pairs$sex1) & !is.na(pairs$sex2)),]

  pairs_unrelated <- subset(pairs_full_sex, pairs_full_sex$kinship<0.05)

  pairs_unrelated_hetero <- pairs_unrelated[-which(pairs_unrelated$sex1==pairs_unrelated$sex2),]
  pairs_unrelated_hetero2 <- pairs_unrelated_hetero[-which(is.na(pairs_unrelated_hetero$sex1) | is.na(pairs_unrelated_hetero$sex2)),]

  temp1 <- pairs_unrelated_hetero[which(pairs_unrelated_hetero$sex1==0),c(2,1,3,4,6,5)]
  colnames(temp1) <- colnames(pairs_unrelated_hetero)

  temp2 <- pairs_unrelated_hetero[which(pairs_unrelated_hetero$sex1==1),]
  pairs_filter <- rbind(temp2, temp1)
  ## 79150      6

  pairs_filter$house_rel_member1 <- household_relationships[[field]][match(pairs_filter[["HOUSEHOLD_MEMBER1"]], household_relationships$userId)]
  pairs_filter$house_rel_member2 <- household_relationships[[field]][match(pairs_filter[["HOUSEHOLD_MEMBER2"]], household_relationships$userId)]
  pairs_filter2 <- pairs_unrelated_hetero[which(pairs_filter$house_rel_member1==TRUE & pairs_filter$house_rel_member2==TRUE ) ,]
  ## 75048     6

  colnames(pairs_filter2)[5] <- "HOUSEHOLD_MEMBER1_sex"
  colnames(pairs_filter2)[6] <- "HOUSEHOLD_MEMBER2_sex"
  return(pairs_filter2)

}

munge_sqc <- function(sqc,fam){
  sqc_sub <- sqc[,c(26:65)]
  colnames(sqc_sub) <- c(sapply(1:40, function(x) {paste0("PC_",x)}))
  sqc_sub$ID <- fam[,1]
  return(sqc_sub)
}

add_pcs <- function(pairs,pheno,sqc_sub){
  out <- list()
  for(i in 1:dim(pairs)[1])
  {
    if(pairs[["HOUSEHOLD_MEMBER1_sex"]][i]==0)
    {
      male_ID <-pairs[["HOUSEHOLD_MEMBER2"]] [i]
      female_ID <-pairs[["HOUSEHOLD_MEMBER1"]] [i]
      pairs[["HOUSEHOLD_MEMBER1"]][i] <-male_ID
      pairs[["HOUSEHOLD_MEMBER2"]][i] <-female_ID
    }

  }
  pairs[["HOUSEHOLD_MEMBER1_sex"]] <- 1
  pairs[["HOUSEHOLD_MEMBER2_sex"]] <- 0

  for(merge_by in c("HOUSEHOLD_MEMBER1", "HOUSEHOLD_MEMBER2"))
  {

    merge1_temp <- merge(pairs, pheno, by.x=merge_by, by.y="userId")
    #length(which(!pairs[,1] %in% pheno1_sub[,1]))
    ### when merge_by="HOUSEHOLD_MEMBER1" -> 3 people were in household file and not pheno
    #length(which(!pairs[,2] %in% pheno1_sub[,1]))
    ### when merge_by="HOUSEHOLD_MEMBER2" -> 1 person were in household file and not pheno

    merge2_temp <- merge(merge1_temp, sqc_sub, by.x=merge_by, by.y="ID")
    #length(which(!merge1_temp[,"HOUSEHOLD_MEMBER1"] %in% sqc_sub[,"ID"]))
    ### when merge_by="INDEX" -> 2009 people were in household file and not geno PC file
    #length(which(!merge1_temp[,"HOUSEHOLD_MEMBER2"] %in% sqc_sub[,"ID"]))
    ### when merge_by="HOUSEHOLD_MEMBER" -> 2799 people were in household file and not geno PC file

    colnames(merge2_temp) <- c(colnames(merge2_temp)[1:6], paste0(merge_by, "_", colnames(merge2_temp)[-c(1:6)]))

    out[[paste0(merge_by, "_pheno_data")]] <- merge2_temp
    #write.csv(merge2_temp, paste0("pipeline/data_setup/", merge_by,"_pheno_model_adjustments", ".csv"),row.names=F, quote=T )

  }

return(out)
}


calc_time_together <- function(pheno_cov,time_at_household,time_at_household_raw){

  households_temp <- pheno_cov[,1:3]

  households_temp$house_time_member1 <- time_at_household[[time_at_address_field]][match(households_temp[["HOUSEHOLD_MEMBER1"]], time_at_household$userId)]
  households_temp$house_time_member2 <- time_at_household[[time_at_address_field]][match(households_temp[["HOUSEHOLD_MEMBER2"]], time_at_household$userId)]
  households_temp$house_time_raw_member1 <- time_at_household_raw[[time_at_address_raw_field]][match(households_temp[["HOUSEHOLD_MEMBER1"]], time_at_household_raw$eid)]
  households_temp$house_time_raw_member2 <- time_at_household_raw[[time_at_address_raw_field]][match(households_temp[["HOUSEHOLD_MEMBER2"]], time_at_household_raw$eid)]

  households_temp <- transform(households_temp, time_together = pmin(house_time_member1, house_time_member2))
  households_temp <- transform(households_temp, time_together_raw = pmin(house_time_raw_member1, house_time_raw_member2))
  return(households_temp)
}


### File A1 ----
## Description: tests phenotypic correlations between household pairs that have been
## filtered for kinship < 0.05, and to be opposite-sex pairs
## run within project folder: /data/sgg2/jenny/projects/MR_Shared_Environment

compute_trait_corr <- function(phesant_directory,UKBB_directory,pairs_filter){
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

SGG_link_with_Neale <- function(Neale_SGG_dir){

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
  return(Neale_SGG_dir_filt2)
}

filter_by_corr <- function(trait_corr,Neale_SGG_dir_filt2,household_correlation_threshold){

  merge_temp <- merge(trait_corr, Neale_SGG_dir_filt2, by.x="ID", by.y="sgg_phesant_name", fill=T)
  merge_temp$r2 <- as.numeric(as.character(merge_temp$r2))
  corr_traits <- merge_temp[which(sqrt(merge_temp$r2) > household_correlation_threshold),]
  corr_traits <- corr_traits[,-which(colnames(corr_traits) %in% c("description.x"))]
  colnames(corr_traits) <- c("SGG_PHESANT_ID", "SGG_PHESANT_ID_sub","N_pairs","r2",
              "Neale_pheno_ID","Neale_pheno_ID_sub", "Neale_file_sex","description",
              "variable_type", "SGG_request", "PHESANT_processed", "SGG_location",
              "Neale_downloaded","category","define_category:T/F",
              "v2_exists", "v2_downloaded")
  return(corr_traits)
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

  return(list("download_rest"=download_rest, "define_cats"=define_cats2))
}

### File A3 ----

## Description: download missing neale files using categories defined here:
## "output/tables/define_Neale_categories_filled.csv".
## Update categories in phenotypic correlations file: output/tables/pheno_household_filter_corr.csv
## run within project folder: /data/sgg2/jenny/projects/MR_Shared_Environment

download_Neale <- function(filled_cats,download_rest,traits_corr_filter, reference_file){
  full_dl_list <- rbind(filled_cats, download_rest)
  select <- dplyr::select

  for(category
     in c("body", "brain", "diet", "disease", "disease_proxy", "lifestyle", "parental_pheno"))
  {
    IDs <- as.character(full_dl_list[which(full_dl_list[,"category"]==category),"Neale_pheno_ID"])
    download_neale_files( IDs,
                      category = category , reference_file=reference_file)
  }

  corr_traits <- traits_corr_filter
  levels(corr_traits$category) <- union (levels(corr_traits$category), levels(filled_cats$category))
  index <- which(corr_traits[["define_category:T/F"]]==TRUE)
  pheno_missing_cat <- corr_traits[index,"Neale_pheno_ID"]
  corr_traits[index, "category"] <- as.factor(filled_cats[["category"]][match(pheno_missing_cat, filled_cats$Neale_pheno_ID)])

  cat(paste0("Missing Neale files have been successfully downloaded and categorized,\n",
  "now they need to be processed (clumped for LD and extracted for p<0.1).\n"))

  return(corr_traits)
}

update_download_info <- function(corr_traits,Neale_SGG_dir){

  corr_traits_both <- corr_traits[which(corr_traits[["Neale_file_sex"]]=="both"),]
  irnt=TRUE

  index <- which(corr_traits[["v2_exists"]]==TRUE)
  index_df <- corr_traits[index,c("Neale_pheno_ID","Neale_file_sex")]
  index_df2 <- merge(index_df, Neale_SGG_dir, by.x=c("Neale_pheno_ID","Neale_file_sex"), by.y=c("phenotype","phenotype_sex"),all.x=TRUE)
  now_downloaded <- index_df2$v2_downloaded
  corr_traits[index, "v2_downloaded"] <- now_downloaded
  return(corr_traits)

}



get_IV_list <- function(corr_traits, i, reference_file) #
{
  errorFile <- paste0("get_IV_list", "_errors.txt")

  corr_traits <- read.csv("output/tables/2.household_correlations.corr_filter.csv", header=T, check.names=F)
  corr_traits_both <- corr_traits[which(corr_traits[["Neale_file_sex"]]=="both"),]


  existing_files_full <- list()
  existing_files_full  =  list.files( Neale_output_path,
                                     recursive = TRUE,
                                     pattern = '[.]gz' )
  IVs_full <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ), full.names=T)
  IVs <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ))
  irnt=TRUE

  tryCatch(
  {


    category <- corr_traits_both[i,"category"]
    Neale_id <- corr_traits_both[i,"Neale_pheno_ID"]
    ID <- corr_traits_both[i,"Neale_pheno_ID"]
    phenotype_ids  =  paste0( '^', ID, ifelse( irnt, '(_irnt|)$', '(_raw|)$' ) ) %>%
        paste( collapse = '|' )

    file_name_temp  =  reference_file %>%
        filter( str_detect( .$'Phenotype Code', phenotype_ids ) & Sex=="both_sexes")

    file_name <- file_name_temp[["Phenotype Code"]][1]
    folder <- IVs_full[which(grepl(paste0("/",file_name ,"\\."), IVs_full) & grepl("both_sexes", IVs))]
    #folder <- IVs_full[which(grepl(paste0(file_name ), IVs_full) & grepl("both_sexes", IVs))]
    if(length(folder)==2){folder <- folder[grep("v2",folder)]} ## use version 2 if it exists
    name <- IVs[which(grepl(paste0("^", file_name ,"\\."), IVs) & grepl("both_sexes", IVs))]
    if(length(name)==2){name <- name[grep("v2",name)]} ## use version 2 if it exists

    IV_snp_list <- numeric()
    for(threshold in c(IV_threshold,GRS_thresholds))
    {
    if(file.exists(paste0( project_dir,"/analysis/data_setup/IV_lists/", Neale_id, "_IVs_", threshold,"_both_sexes.txt")))
      {file.remove(paste0( project_dir,"/analysis/data_setup/IV_lists/", Neale_id, "_IVs_", threshold,"_both_sexes.txt"))}
    }
    threshold <- IV_threshold ## only extract IVs
    for(chr in 1:22)

    {

        original_data <- read.table(paste0(folder,"/", name,"_unpruned_chr", chr, ".txt"), header=T)
        if(file.exists(paste0(folder,"/", name,"_unpruned_chr", chr, ".clumped")))
        {
          clump_data <- read.table(paste0(folder,"/", name,"_unpruned_chr", chr, ".clumped"),he=T)
          clumped_SNPs <- clump_data[,"SNP"]
          IV_temp <- original_data[which(original_data[["P"]]<threshold & original_data[["SNP"]] %in% clumped_SNPs),"SNP"]
          write.table(IV_temp, paste0( project_dir,"/pipeline/data_setup/IV_lists/", Neale_id, "_IVs_", threshold,"_both_sexes.txt"), append=T, row.names=F, col.names=F, quote=F)
        }

    }


  },
    error = function(e) {cat(c(i, paste0(as.character(e),"\n")), file=paste0(project_dir,"/",errorFile), append=T)},
    warning = function(w) {cat(c(i, paste0(as.character(w),"\n")), file=paste0(project_dir,"/",errorFile), append=T)}
  )

}



summarize_IVs <- function(corr_traits, i, reference_file)
{

  IVs_full <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ), full.names=T)
  IVs <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ))
  category <- corr_traits_both[i,"category"]
  Neale_id <- corr_traits_both[i,"Neale_pheno_ID"]
  ID <- corr_traits_both[i,"Neale_pheno_ID"]
  phenotype_ids  =  paste0( '^', ID, ifelse( irnt, '(_irnt|)$', '(_raw|)$' ) ) %>%
      paste( collapse = '|' )

  file_name_temp  =  reference_file %>%
      filter( str_detect( .$'Phenotype Code', phenotype_ids ) & Sex=="both_sexes")

  file_name <- file_name_temp[["Phenotype Code"]]
  folder <- IVs_full[which(grepl(file_name, IVs_full) & grepl("both_sexes", IVs))]
  name <- IVs[which(grepl(file_name, IVs) & grepl("both_sexes", IVs))]
  IV_snp_list <- numeric()
  threshold <- IV_threshold
  #if(file.exists(paste0( "pipeline/data_setup/IV_lists/", Neale_id, "_IVs_", threshold,"_both_sexes.txt")))
  if(file.info(paste0( "analysis/data_setup/IV_lists/", Neale_id, "_IVs_", threshold,"_both_sexes.txt"))$size!=0)
  {
    GW_ivs <- fread( paste0( "analysis/data_setup/IV_lists/", Neale_id, "_IVs_", threshold,"_both_sexes.txt"),data.table=F, header=F)
    out <- cbind(i,dim(GW_ivs)[1])
  } else out <- cbind(i,0)

  return(out)

}

IV_filter <- function(corr_traits, IV_summary, num_IVs_threshold){

  corr_traits_both <- corr_traits[which(corr_traits[["Neale_file_sex"]]=="both"),]
  order_summary <- IV_summary[order(IV_summary[,1]),]
  corr_traits_both$num_IVs <- order_summary[,2]

  to_run <- corr_traits_both[which(corr_traits_both$num_IVs >=num_IVs_threshold),]
  return(list(to_run = to_run, non_filtered = corr_traits_both))

  cat(paste0("Household phenotypic correlations have been filtered to include only phenotypes
  which have at least minimum number of IVs according to specified threshold
  (see 'scripts/settings.r' for chosen threshold),\n",
  "filtered phenotypic correlations have been saved to:
  'output/tables/3.household_correlations.numIVs_filter.csv'.\n"))

}

calc_sex_het <- function(traits,i,variant_data,output_folder){

  irnt=TRUE
  existing_files_full <- list()
  existing_files_full  =  list.files( Neale_output_path,
                                     recursive = TRUE,
                                     pattern = '[.]gz' )


  IVs_full <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ), full.names=T)
  IVs <- list.files(path=paste0(Neale_summary_dir,"/IVs/clump/" ))

  #traits$num_IVs_pass_het <- NA

  result <- NA

  trait_ID <- as.character(traits[i,"Neale_pheno_ID"])
  phenotype_ids  =  paste0( '^', trait_ID, ifelse( irnt, '(_irnt|)$', '(_raw|)$' ) ) %>%
     paste( collapse = '|' )
  file_name_temp  =  reference_file %>%
     filter( str_detect( .$'Phenotype Code', phenotype_ids ) & Sex=="both_sexes")
  file_name <- file_name_temp[["Phenotype Code"]][1]

  #trait_info <- read.table(paste0(pheno_dir,"/trait_info.txt"), header=F, row.names=1)

  IV_list_both_sexes <- fread(paste0( project_dir,"/pipeline/data_setup/IV_lists/", trait_ID, "_IVs_", IV_threshold,"_both_sexes.txt"), data.table=F, header=F)
  SNP_rows <- which(variant_data[,"rsid"] %in% IV_list_both_sexes[,1])
  IV_folder <- IVs_full[which(grepl(paste0("/",file_name ,"\\."), IVs_full) & grepl("both_sexes", IVs))]
  IV_file_name <- IVs[which(grepl(paste0("^", file_name ,"\\."), IVs) & grepl("both_sexes", IVs))]

  if(length(IV_folder)==2){IV_folder <- IV_folder[grep("v2",IV_folder)]} ## use version 2 if it exists
  if(length(IV_file_name)==2){IV_file_name <- IV_file_name[grep("v2",IV_file_name)]} ## use version 2 if it exists

  for(exposure_sex in c("both_sexes", "male", "female"))
  {

    specific_IV_folder <- IVs_full[which(grepl(paste0("/",file_name ,"\\."), IVs_full) & grepl(paste0("[.]",exposure_sex), IVs))]
    specific_IV_name <- IVs[which(grepl(paste0("^", file_name ,"\\."), IVs) & grepl(paste0("[.]",exposure_sex), IVs))]
    if(length(specific_IV_folder)==2){specific_IV_folder <- specific_IV_folder[grep("v2",specific_IV_folder)]} ## use version 2 if it exists
    if(length(specific_IV_name)==2){specific_IV_name <- specific_IV_name[grep("v2",specific_IV_name)]} ## use version 2 if it exists

    original_neale_temp <- existing_files_full[which(grepl(paste0("/",gsub(".IVs", "",specific_IV_name)), existing_files_full))]
    original_neale_file <- paste0(Neale_output_path, "/", original_neale_temp)
    assign(paste0(exposure_sex, "_IV_folder"), specific_IV_folder)
    assign(paste0(exposure_sex, "_original_Neale_file"), original_neale_file)

  }


  for(file in c("male", "female"))
  {

    temp <- fread(get(paste0(file,"_original_Neale_file")), data.table=F)
    temp <- temp[SNP_rows,]
    exp_var_merge <- cbind(temp,variant_data[SNP_rows,])
    IV_cols <- c("rsid", "chr", "beta", "se", "pval", "ref", "alt", "AF","n_complete_samples")
    exp_var_merge <- exp_var_merge[,IV_cols]

    assign(paste0(file, "_IV_data"), exp_var_merge)

  }

  IV_list_both_sexes$p_het <- NA
  for(snp in 1:length(SNP_rows))
  {
    beta_F <-  female_IV_data[snp,"beta"]
    beta_M <-  male_IV_data[snp,"beta"]
    SE_F <- female_IV_data[snp,"se"]
    SE_M <- male_IV_data[snp,"se"]
    n_M <- male_IV_data[snp,"n_complete_samples"]
    n_F <- female_IV_data[snp,"n_complete_samples"]
    #  se <- sqrt( (SE_F^2/n_F) + (SE_M^2/n_M) )
    se <- sqrt( (SE_F^2) + (SE_M^2) )
    t <- (beta_F-beta_M)/se
    p_het <- 2*pnorm(-abs(t))
    IV_list_both_sexes$p_het[snp] <- p_het
  }



  IV_list_filter <- IV_list_both_sexes[which(IV_list_both_sexes$p_het > 0.05/length(SNP_rows)),]

  num_pass_filter <- length(which(IV_list_both_sexes$p_het > 0.05/length(SNP_rows)))


  colnames(IV_list_both_sexes) <- c("SNP", "P-het")
  write.table(IV_list_both_sexes, paste0( output_folder, trait_ID, "_sex_het.txt"), row.names=F, col.names=T, quote=F)

  char_row <- data.frame(lapply(traits[i,], as.character), stringsAsFactors=FALSE)

  return(unlist(c(char_row, num_pass_filter)))

}

sex_het_filter <- function(corr_traits, sex_het_summary, num_IVs_threshold){

  colnames(sex_het_summary) <- c(colnames(corr_traits),"num_IVs_pass_het")
  result_df <- as.data.frame(sex_het_summary)
  result_df$num_IVs_pass_het <- as.numeric(as.character(result_df$num_IVs_pass_het))
  to_run <- result_df[which(result_df$num_IVs_pass_het >=num_IVs_threshold),]
  removed <- result_df[-which(result_df$num_IVs_pass_het >=num_IVs_threshold),]
  return(list(to_run = to_run, non_filtered = result_df))

  cat(paste0("Household phenotypic correlations have been filtered to include only phenotypes
  which have at least minimum number of IVs according to specified threshold after removing IVs with
  with significant evidence of heterogeneity between sexes (p<0.05/[number of GW-sig IVs]),
  filtered phenotypic correlations have been saved to:
  'output/tables/4.household_correlations.sexhet_filter.csv'.\n"))
}
