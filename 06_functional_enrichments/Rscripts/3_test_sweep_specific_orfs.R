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

sister_differentiating_genes <- read_csv("intermediate_files/genes_differentiating_sweeps_from_sisters.csv") %>%
  pull(qseqid)

sweep_specific_orfs_all_df <- read_csv("intermediate_files/unique_prots_across_sweep.csv") %>%
  # filter for genes unique to sister sweeps
  filter(ORF_id %in% sister_differentiating_genes) %>%
  filter(COG_category != "-") %>%
  filter(COG_category != "S") %>%
  
  # only test single categories
  filter(nchar(COG_category) == 1)

write.csv(sweep_specific_orfs_all_df, "intermediate_files/final_sweep_specific_ORFs.csv", row.names = FALSE)

# loop through all SGBs
SGB_vector <- read_csv("meta_data/SGB_list.csv") %>%
  pull(SGB)

j <- 1
list_collect_cog_pvalues <- list()

for (SGB in SGB_vector){
  
  # filter for genes unique to sister sweeps
  sweep_specific_orfs_df <- sweep_specific_orfs_all_df %>%
    filter(grepl(SGB, Clonalframe))
  
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
  
  for (i in 1:nrow(cog_count_combined_df)) {
    
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
    
    p_values[i] <- test$p.value
  }
  
  results_cog <- data.frame(COG_category = cog_count_combined_df$COG_category, p_value = p_values)
  results_cog$adj_p_values <- p.adjust(p_values, method = "bonferroni")
  
  results_cog_select <- results_cog %>%
    filter(adj_p_values < 0.000005)
  
  # prepare data for pfam test
  
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
    group_by(PFAM) %>%
    mutate(count_ss = n()) %>%
    ungroup() %>%
    distinct()
  
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
  
  for (i in 1:nrow(pfam_count_combined_df)) {
    
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
    
    p_values[i] <- test$p.value
  }
  
  results_pfam <- data.frame(PFAM = pfam_count_combined_df$PFAM, p_value = p_values)
  results_pfam$adj_p_values <- p.adjust(p_values, method = "bonferroni")
  
  results_pfam_select <- results_pfam %>%
    filter(adj_p_values < 0.000005)
  
  #
  
  
  # write the enriched COG categories and PFAMs to file
  if(nrow(results_cog_select)){
    results_cog_select$SGB <- SGB
    write_cogs <- paste("sig_sweep_specific_orfs/", SGB, "_enriched_COG_categories_000005.csv", sep = "")
    write.csv(results_cog_select, write_cogs, row.names = FALSE)
  }
  
  if(nrow(results_pfam_select)){
    results_pfam_select$SGB <- SGB
    write_pfams <- paste("sig_sweep_specific_orfs/", SGB, "_enriched_PFAMs_000005.csv", sep = "")
    write.csv(results_pfam_select, write_pfams, row.names = FALSE)
  }
  
  results_cog$SGB <- SGB
  list_collect_cog_pvalues[[j]] <- results_cog
  j <- j + 1
}

df_save_cog_pvalues <- bind_rows(list_collect_cog_pvalues)
write.csv(df_save_cog_pvalues, "sig_sweep_specific_orfs/pvalues_COG_unique_genes.csv")