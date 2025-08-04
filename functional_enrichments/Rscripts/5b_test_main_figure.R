library(readr)
library(tidyr)
library(dplyr)


proka_clonalframes_df <- read_delim("raw_data/clonalframes.tsv", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)

clonalframe_map <- proka_clonalframes_df %>%
  select(locus_tag, clonalframe) %>%
  rename(ORF_id = locus_tag) %>%
  mutate(clonalframe_clean = gsub("_cf_consensus", "", clonalframe)) %>%
  select(-clonalframe)

clonalframe_annotations_all_df <- read_delim("raw_data/clonalframes.emapper.annotations", 
                                             delim = "\t", escape_double = FALSE, 
                                             trim_ws = TRUE, skip = 4)
colnames(clonalframe_annotations_all_df)[1] <- "ORF_id"

clonalframe_annotations_all_df <- clonalframe_annotations_all_df %>%
  inner_join(clonalframe_map)

######

sweep_specific_orfs_all_df <- read_csv("intermediate_files/final_sweep_specific_ORFs.csv") %>%
  mutate(Sweep = gsub("_cf_consensus", "", Clonalframe)) %>%
  mutate(Sweep = gsub("diffH", "", Sweep)) 

main_figure_COG <- read_csv("intermediate_files/main_figure_COG.csv") %>%
  rename(Annotation_category = COG_category)

df_disease_associations <- read_csv("meta_data/associated_sweeps_unique_genes.csv") %>%
  select(Sweep, Enrichment) %>% 
  mutate(Enrichment = if_else(grepl("enriched", Enrichment), "Enriched", 
                                                                   if_else(grepl("depleted", Enrichment), "Depleted", "Other"))) %>%
  distinct() 

sweep_specific_orfs_all_df <- left_join(sweep_specific_orfs_all_df, df_disease_associations)

sweep_specific_orfs_all_df$Enrichment[is.na(sweep_specific_orfs_all_df$Enrichment)] <- "None"

SGB_vector <- main_figure_COG %>%
  pull(SGB)

associations_vector <- df_disease_associations %>%
  pull(Enrichment)

associations_vector <- c(associations_vector, "None")

j <- 1
list_collect_cog_pvalues <- list()

df_output <- data.frame(SGB_col = character(0), Enrichment_col = character(0), Pvalue_col = numeric(0), stringsAsFactors = FALSE)

for (SGB in SGB_vector){
  for (Assoc in associations_vector){

    # filter for genes unique to sister sweeps
    sweep_specific_orfs_df <- sweep_specific_orfs_all_df %>%
      filter(grepl(SGB, Clonalframe)) %>%
      filter(grepl(Assoc, Enrichment))
    
    if(nrow(sweep_specific_orfs_df) == 0){ next }
    
    # prepare data for cog test
    
    cog_counts_ss_df <- sweep_specific_orfs_df %>% 
      filter(!(is.na(COG_category))) %>%
      select(COG_category) %>%
      group_by(COG_category) %>%
      mutate(count_ss = n()) %>%
      ungroup() %>%
      distinct() %>%
      filter(COG_category != "-")
    
    clonalframe_annotations_df <- clonalframe_annotations_all_df %>%
      filter(grepl(SGB, clonalframe_clean))
    
    cog_counts_cf_df <- clonalframe_annotations_df %>%
      separate(ORF_id, into = c("prokka_id"), sep = "_", remove = FALSE) %>% 
      filter(!(is.na(COG_category))) %>%
      select(COG_category) %>%
      group_by(COG_category) %>%
      mutate(count_cf = n()) %>%
      ungroup() %>%
      distinct() %>%
      filter(COG_category != "-") %>%
      filter(COG_category != "S") %>%
      
      # only test single categories
      filter(nchar(COG_category) == 1)
    
    cog_count_combined_df <- left_join(cog_counts_cf_df, cog_counts_ss_df)
    cog_count_combined_df$count_ss[is.na(cog_count_combined_df$count_ss)] <- 0
    
    # test COG categories
    
    total_count_cf <- sum(cog_count_combined_df$count_cf)
    total_count_ss <- sum(cog_count_combined_df$count_ss)
    
    p_values <- numeric(length(cog_count_combined_df$COG_category))
    
    i <- which(cog_count_combined_df$COG_category == "M")
    
    cog_category_count_cf <- cog_count_combined_df$count_cf[i]
    cog_category_ss <- cog_count_combined_df$count_ss[i]
    
    other_categories_count_cf <- total_count_cf - cog_category_count_cf
    other_categories_count_ss <- total_count_ss - cog_category_ss
    
    cog_category_count_not_ss = cog_category_count_cf - cog_category_ss
    other_categories_count_not_ss = other_categories_count_cf - other_categories_count_ss
    
    dat <- data.frame(
      "in_sweep_specific" = c(cog_category_ss, other_categories_count_ss),
      "not_in_sweep_specific" = c(cog_category_count_not_ss, other_categories_count_not_ss),
      row.names = c("Cog-Category", "All-Other_Categories"),
      stringsAsFactors = FALSE
    )
    colnames(dat) <- c("In-Sweep-Specific-ORFs", "Not-In-Sweep-Specific-ORFs")
    
    test <- fisher.test(dat, alternative = "greater")
      
    p_value <- test$p.value
    
    new_row <- data.frame(SGB_col = SGB, Enrichment = Assoc, Pvalue_col = p_value,stringsAsFactors = FALSE)
    df_output <- rbind(df_output, new_row) %>%
      distinct()
    
    }
}

write.csv(df_output, "sig_sweep_specific_orfs/pvalues_main_figure_COG.csv", row.names = FALSE)

######
j <- 1
list_collect_cog_pvalues <- list()

df_output <- data.frame(SGB_col = character(0), Enrichment_col = character(0), Pvalue_col = numeric(0), stringsAsFactors = FALSE)

for (SGB in SGB_vector){
  for (Assoc in associations_vector){
    
    # filter for genes unique to sister sweeps
    sweep_specific_orfs_df <- sweep_specific_orfs_all_df %>%
      filter(grepl(SGB, Clonalframe)) %>%
      filter(grepl(Assoc, Enrichment))
    
    if(nrow(sweep_specific_orfs_df) == 0){ next }

    pfam_counts_ss_df <- sweep_specific_orfs_df %>%
      filter(!(is.na(PFAMs)), PFAMs != "-") %>%
      select(PFAMs) %>%
      mutate(row_id = row_number()) %>%
      separate_rows(PFAMs, sep = ",") %>%
      group_by(row_id) %>%
      mutate(col_name = paste0("PFAM_", row_number())) %>%
      ungroup() %>%
      pivot_wider(names_from = col_name, values_from = PFAMs) %>%
      select(-row_id) %>%
      gather(colnames, PFAM) %>%
      filter(!(is.na(PFAM))) %>%
      select(-colnames) %>%
      
      mutate(PFAM = if_else(grepl("Glycos_trans|Glyco_trans", PFAM), "Glyco_trans", 
                          if_else(grepl("DDE_Tnp", PFAM), "DDE_Tnp", 
                          if_else(grepl("Phage_integrase", PFAM), "Phage_integrase", 
                          if_else(grepl("Methylase_S", PFAM), "Methylase_S", 
                          if_else(grepl("HTH_32", PFAM), "HTH_32", 
                          if_else(grepl("rve", PFAM), "rve", 
                          if_else(grepl("RVT_1", PFAM), "RVT_1", "other")))))))) %>%
      
      group_by(PFAM) %>%
      mutate(count_ss = n()) %>%
      ungroup() %>%
      distinct()
    
    clonalframe_annotations_df <- clonalframe_annotations_all_df %>%
      filter(grepl(SGB, clonalframe_clean))
    
    pfam_counts_cf_df <- clonalframe_annotations_df %>%
      separate(ORF_id, into = c("prokka_id"), sep = "_", remove = FALSE) %>%
      filter(!(is.na(PFAMs)), PFAMs != "-") %>%
      select(PFAMs) %>%
      mutate(row_id = row_number()) %>%
      separate_rows(PFAMs, sep = ",") %>%
      group_by(row_id) %>%
      mutate(col_name = paste0("PFAM_", row_number())) %>%
      ungroup() %>%
      pivot_wider(names_from = col_name, values_from = PFAMs) %>%
      select(-row_id) %>%
      gather(colnames, PFAM) %>%
      filter(!(is.na(PFAM))) %>%
      select(-colnames) %>%
      
      mutate(PFAM = if_else(grepl("Glycos_trans|Glyco_trans", PFAM), "Glyco_trans", 
                    if_else(grepl("DDE_Tnp", PFAM), "DDE_Tnp", 
                    if_else(grepl("Phage_integrase", PFAM), "Phage_integrase", 
                    if_else(grepl("Methylase_S", PFAM), "Methylase_S", 
                    if_else(grepl("HTH_32", PFAM), "HTH_32", 
                    if_else(grepl("rve", PFAM), "rve", 
                    if_else(grepl("RVT_1", PFAM), "RVT_1", "other")))))))) %>%
      
      group_by(PFAM) %>%
      mutate(count_cf = n()) %>%
      ungroup() %>%
      distinct()
    
    pfam_count_combined_df <- left_join(pfam_counts_cf_df, pfam_counts_ss_df)
    pfam_count_combined_df$count_ss[is.na(pfam_count_combined_df$count_ss)] <- 0
    
    # testing pfams
    
    total_count_cf <- sum(pfam_count_combined_df$count_cf)
    total_count_ss <- sum(pfam_count_combined_df$count_ss)
    
    p_values <- numeric(length(pfam_count_combined_df$PFAM))
    
    i <- which(pfam_count_combined_df$PFAM == "Glyco_trans")
    
    pfam_category_count_cf <- pfam_count_combined_df$count_cf[i]
    pfam_category_ss <- pfam_count_combined_df$count_ss[i]
    
    other_categories_count_cf <- total_count_cf - pfam_category_count_cf
    other_categories_count_ss <- total_count_ss - pfam_category_ss
    
    pfam_category_count_not_ss = pfam_category_count_cf - pfam_category_ss
    other_categories_count_not_ss = other_categories_count_cf - other_categories_count_ss
    
    dat <- data.frame(
      "in_sweep_specific" = c(pfam_category_ss, other_categories_count_ss),
      "not_in_sweep_specific" = c(pfam_category_count_not_ss, other_categories_count_not_ss),
      row.names = c("PFAM-Category", "All-Other_Categories"),
      stringsAsFactors = FALSE
    )
    colnames(dat) <- c("In-Sweep-Specific-ORFs", "Not-In-Sweep-Specific-ORFs")
    
    test <- fisher.test(dat, alternative = "greater")
    
    p_value <- test$p.value
  
    new_row <- data.frame(SGB_col = SGB, Enrichment = Assoc, Pvalue_col = p_value,stringsAsFactors = FALSE)
    df_output <- rbind(df_output, new_row) %>%
      distinct()

  }
}

write.csv(df_output, "sig_sweep_specific_orfs/pvalues_main_figure_PFAM.csv", row.names = FALSE)
