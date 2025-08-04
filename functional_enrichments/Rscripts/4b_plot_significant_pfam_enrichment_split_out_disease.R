library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

input_dir <- "sig_sweep_specific_orfs/"

input_files <- list.files(input_dir, pattern = "PFAMs_000005.csv")

df_tmp_list <- list()

i <- 1

for (file in input_files){
  input_path <- paste(input_dir, file, sep = "")
  df_tmp <- read_csv(input_path, col_types = cols(.default = col_character()))
  df_tmp_list[[i]] <- df_tmp
  i <- i + 1
}

df_test_compile <- bind_rows(df_tmp_list) %>%
  mutate(SGB = gsub("diffH", "", SGB))

#####

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

#####

# Clean up pfams

Phage_integrase <- c("Phage_integrase")
Methylase_S <- c("Methylase_S")
HTH_32 <- c("HTH_32")
rve <- c("rve")
RVT_1 <- c("RVT_1")

#####

df_percent_all_genes <- clonalframe_annotations_all_df %>%
  
  # clean up pfams
  mutate(PFAM_clean = if_else(grepl("Glycos_trans|Glyco_trans", PFAMs), "Glyco_trans", 
                      if_else(grepl("DDE_Tnp", PFAMs), "DDE_Tnp", 
                      if_else(grepl("Phage_integrase", PFAMs), "Phage_integrase", 
                      if_else(grepl("Methylase_S", PFAMs), "Methylase_S", 
                      if_else(grepl("HTH_32", PFAMs), "HTH_32", 
                      if_else(grepl("rve", PFAMs), "rve", 
                      if_else(grepl("RVT_1", PFAMs), "RVT_1", "other")))))))) %>%
  
  select(clonalframe_clean, ORF_id, PFAM_clean) %>%
  
  separate(clonalframe_clean, into = c("SGB"), sep = "_", remove = FALSE) %>%
  
  mutate(SGB = gsub("diffH", "", SGB)) %>%
  
  group_by(SGB) %>%
  mutate(Total = n()) %>%
  ungroup() %>%
  
  group_by(SGB, PFAM_clean) %>%
  mutate(Total_per_PFAM = n()) %>%
  ungroup() %>%
  
  select(SGB, PFAM_clean, Total, Total_per_PFAM) %>%
  
  mutate(Per_all_PFAM = (Total_per_PFAM / Total)*100) %>%
  distinct() %>%
  
  select(SGB, PFAM_clean, Per_all_PFAM) %>%
  
  filter(PFAM_clean != "other")

####

#####

df_disease_associations <- read_csv("meta_data/associated_sweeps_unique_genes.csv") %>%
  select(Sweep, Enrichment) %>%
  distinct()

df_final_sweep_specific_ORFs <- read_csv("intermediate_files/final_sweep_specific_ORFs.csv") %>%
  mutate(Sweep = gsub("_cf_consensus", "", Clonalframe)) %>%
  mutate(Sweep = gsub("diffH", "", Sweep)) 

df_percent_sweep_specific_genes <- df_final_sweep_specific_ORFs %>%
  
  mutate(PFAM_clean = if_else(grepl("Glycos_trans|Glyco_trans", PFAMs), "Glyco_trans", 
                      if_else(grepl("DDE_Tnp", PFAMs), "DDE_Tnp", 
                      if_else(grepl("Phage_integrase", PFAMs), "Phage_integrase", 
                      if_else(grepl("Methylase_S", PFAMs), "Methylase_S", 
                      if_else(grepl("HTH_32", PFAMs), "HTH_32", 
                      if_else(grepl("rve", PFAMs), "rve", 
                      if_else(grepl("RVT_1", PFAMs), "RVT_1", "other")))))))) %>%
  
  separate(Clonalframe, into = c("SGB"), sep = "_", remove = FALSE) %>%
  
  mutate(SGB = gsub("diffH", "", SGB)) %>%
  
  left_join(df_disease_associations) %>%
  
  mutate(Pos_assoc_disease = if_else(grepl("enriched", Enrichment), "Enriched", 
                                     if_else(grepl("depleted", Enrichment), "Depleted", "Other"))) %>%
  
  
  select(SGB, Clonalframe, PFAM_clean, Pos_assoc_disease, Sweep) %>% 
  
  
  group_by(SGB) %>%
  mutate(Total_sweep_specific = n()) %>%
  ungroup() %>%
  
  group_by(SGB, Pos_assoc_disease) %>%
  mutate(Total_sweep_specific_per_association = n()) %>%
  ungroup() %>%
  
  group_by(SGB, PFAM_clean) %>%
  mutate(Total_sweep_specific_per_PFAM = n()) %>%
  ungroup() %>%
  
  group_by(SGB, PFAM_clean, Pos_assoc_disease) %>%
  mutate(Total_sweep_specific_per_PFAM_per_association = n()) %>%
  ungroup() %>%
  
  select(SGB, Pos_assoc_disease, PFAM_clean, Total_sweep_specific, Total_sweep_specific_per_PFAM, Total_sweep_specific_per_PFAM_per_association, Total_sweep_specific_per_association) %>%
  
  mutate(Per_sweep_specific_per_PFAM = (Total_sweep_specific_per_PFAM / Total_sweep_specific)*100) %>%
  mutate(Per_sweep_specific_per_PFAM_per_association = (Total_sweep_specific_per_PFAM_per_association / Total_sweep_specific_per_association)*100) %>%
  
  distinct() %>%
  
  select(SGB, Pos_assoc_disease, PFAM_clean, Per_sweep_specific_per_PFAM, Per_sweep_specific_per_PFAM_per_association) %>%
  
  filter(PFAM_clean != "other")

###

df_taxa <- read_csv("meta_data/TableS2_sweep_summary_values.csv") %>%
  
  select(SGB, Category, Family, Genus, Species)

df_plot <- inner_join(df_test_compile, df_percent_all_genes) %>%
  inner_join(df_percent_sweep_specific_genes) %>%
  inner_join(df_taxa) %>%
  distinct() %>%
  filter(Category == "commensals") %>%
  
  select(-p_value, -adj_p_values, -PFAM, -Category) %>%
  filter(PFAM_clean == "Glyco_trans" | PFAM_clean == "Phage_integrase" |PFAM_clean == "DDE_Tnp" |PFAM_clean == "Methylase_S")

df_plot_all <- df_plot %>%
  select(-Pos_assoc_disease, -Per_sweep_specific_per_PFAM_per_association) %>%
  distinct() %>%
  
  gather(Group, Percent, -PFAM_clean, -SGB, -Family, -Genus, -Species) 


# ggplot(data = df_plot_all, aes(x = Species, y = Percent, fill = Group)) +
#   
#   facet_wrap(~PFAM_clean, ncol = 8, scales = "free_x") +
#   
#   geom_bar(stat = "identity", position = position_dodge()) +
#   
#   theme_minimal() + 
#   
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#####

df_plot_disease_asso <- df_plot %>%
  
  group_by(SGB) %>%
  mutate(n_cats = length(unique(Pos_assoc_disease))) %>%
  ungroup() %>%
  filter(n_cats >= 2) %>%
  
  
  
  select(-Per_sweep_specific_per_PFAM, -n_cats) %>%
  distinct() %>%
  
  gather(Group, Percent, -Pos_assoc_disease, -PFAM_clean, -SGB, -Family, -Genus, -Species) 

df_plot_disease_asso_v2 <- df_plot_disease_asso %>%
  mutate(Combined_ID = if_else(Group == "Per_all_PFAM", "All_Sweeps", paste("Sweep_specific_", Pos_assoc_disease, sep = ""))) %>%
  select(PFAM_clean, SGB, Species, Percent, Combined_ID) %>%
  distinct()

# ggplot(data = df_plot_disease_asso_v2, aes(x = Species, y = Percent, fill = Combined_ID)) +
#   
#   facet_grid(.~PFAM_clean, scales = "free_x") +
#   
#   geom_bar(stat = "identity", position = position_dodge()) +
#   
#   theme_minimal() + 
#   
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


main_figure_PFAM_table <- df_plot_disease_asso_v2 %>%
  filter(PFAM_clean == "Glyco_trans")

write.csv(main_figure_PFAM_table, "intermediate_files/main_figure_PFAM.csv", row.names = FALSE)
write.csv(df_plot_disease_asso_v2, "intermediate_files/supp_figure_PFAM.csv", row.names = FALSE)
