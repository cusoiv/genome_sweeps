library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

input_dir <- "sig_sweep_specific_orfs/"

input_files <- list.files(input_dir, pattern = "COG_categories_000005.csv")

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

df_percent_all_genes <- clonalframe_annotations_all_df %>%
  
  filter(COG_category != "-") %>%
  filter(COG_category != "S") %>%
  filter(nchar(COG_category) == 1) %>%
  
  select(clonalframe_clean, ORF_id, COG_category) %>%
  
  separate(clonalframe_clean, into = c("SGB"), sep = "_", remove = FALSE) %>%

  mutate(SGB = gsub("diffH", "", SGB)) %>%
  
  group_by(SGB) %>%
  mutate(Total = n()) %>%
  ungroup() %>%
  
  group_by(SGB, COG_category) %>%
  mutate(Total_per_COG = n()) %>%
  ungroup() %>%
  
  select(SGB, COG_category, Total, Total_per_COG) %>%
  
  mutate(Per_all_COG = (Total_per_COG / Total)*100) %>%
  distinct() %>%
  
  select(SGB, COG_category, Per_all_COG)

####

df_disease_associations <- read_csv("meta_data/associated_sweeps_unique_genes.csv") %>%
  select(Sweep, Enrichment) %>%
  distinct()

df_final_sweep_specific_ORFs <- read_csv("intermediate_files/final_sweep_specific_ORFs.csv") %>%
  mutate(Sweep = gsub("_cf_consensus", "", Clonalframe)) %>%
  mutate(Sweep = gsub("diffH", "", Sweep)) 

df_percent_sweep_specific_genes <- df_final_sweep_specific_ORFs %>%
  
  separate(Clonalframe, into = c("SGB"), sep = "_", remove = FALSE) %>%
  
  mutate(SGB = gsub("diffH", "", SGB)) %>%
  
  left_join(df_disease_associations) %>%
  
  mutate(Pos_assoc_disease = if_else(grepl("enriched", Enrichment), "Enriched", 
                             if_else(grepl("depleted", Enrichment), "Depleted", "Other"))) %>%
  
  select(SGB, Clonalframe, COG_category, Pos_assoc_disease, Sweep) %>% 
  
  group_by(SGB) %>%
  mutate(Total_sweep_specific = n()) %>%
  ungroup() %>%
    
  group_by(SGB, Pos_assoc_disease) %>%
  mutate(Total_sweep_specific_per_association = n()) %>%
  ungroup() %>%
  
  group_by(SGB, COG_category) %>%
  mutate(Total_sweep_specific_per_COG = n()) %>%
  ungroup() %>%
  
  group_by(SGB, COG_category, Pos_assoc_disease) %>%
  mutate(Total_sweep_specific_per_COG_per_association = n()) %>%
  ungroup() %>%
  
  select(SGB, COG_category, Pos_assoc_disease, Total_sweep_specific, Total_sweep_specific_per_association, Total_sweep_specific_per_COG, Total_sweep_specific_per_COG_per_association) %>%
  
  mutate(Per_sweep_specific_per_COG = (Total_sweep_specific_per_COG / Total_sweep_specific)*100) %>%
  mutate(Per_sweep_specific_per_COG_per_association = (Total_sweep_specific_per_COG_per_association / Total_sweep_specific_per_association)*100) %>%
  
  distinct() %>%
  
  select(SGB, Pos_assoc_disease, COG_category, Per_sweep_specific_per_COG, Per_sweep_specific_per_COG_per_association) 
  
###

df_taxa <- read_csv("meta_data/TableS2_sweep_summary_values.csv") %>%
  
  select(SGB, Category, Family, Genus, Species)

df_plot <- inner_join(df_test_compile, df_percent_all_genes) %>%
  inner_join(df_percent_sweep_specific_genes) %>%
  inner_join(df_taxa) %>%
  distinct() %>%
  filter(Category == "commensals") %>%
  
  select(-p_value, -adj_p_values, -Category) 


df_plot_all <- df_plot %>%
  select(-Pos_assoc_disease, -Per_sweep_specific_per_COG_per_association) %>%
  distinct() %>%

  gather(Group, Percent, -COG_category, -SGB, -Family, -Genus, -Species) 
  
  
#df_plot_all$Species <- factor(df_plot$Species, levels = c("s__Bacteroides_fragilis", "s__Bacteroides_cellulosilyticus", "s__Parabacteroides_distasonis", "s__Bacteroides_intestinalis", "s__Bacteroides_ovatus", "s__Bacteroides_uniformis", "s__Phocaeicola_vulgatus", "s__Enterococcus_faecalis"))
#df_plot$COG_category <- factor(df_plot$COG_category, levels = c("M", "L"))


# ggplot(data = df_plot_all, aes(x = Species, y = Percent, fill = Group)) +
#   
#   facet_wrap(~COG_category, ncol = 8, scales = "free_x") +
#   
#   geom_bar(stat = "identity", position = position_dodge()) +
#   
#   theme_minimal() + 
#   
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
######

df_plot_disease_asso <- df_plot %>%
  
  group_by(SGB) %>%
  mutate(n_cats = length(unique(Pos_assoc_disease))) %>%
  ungroup() %>%
  filter(n_cats >= 2) %>%
  
  
  
  select(-Per_sweep_specific_per_COG, -n_cats) %>%
  distinct() %>%
  
  gather(Group, Percent, -Pos_assoc_disease, -COG_category, -SGB, -Family, -Genus, -Species) 


#df_plot_disease_asso$Species <- factor(df_plot$Species, levels = c("s__Bacteroides_fragilis", "s__Bacteroides_cellulosilyticus", "s__Parabacteroides_distasonis", "s__Bacteroides_intestinalis", "s__Bacteroides_ovatus", "s__Bacteroides_uniformis", "s__Phocaeicola_vulgatus", "s__Enterococcus_faecalis"))
#df_plot$COG_category <- factor(df_plot$COG_category, levels = c("M", "L"))

df_plot_disease_asso_v2 <- df_plot_disease_asso %>%
  mutate(Combined_ID = if_else(Group == "Per_all_COG", "All_Sweeps", paste("Sweep_specific_", Pos_assoc_disease, sep = ""))) %>%
  select(COG_category, SGB, Species, Percent, Combined_ID) %>%
  distinct()

# ggplot(data = df_plot_disease_asso_v2, aes(x = Species, y = Percent, fill = Combined_ID)) +
#   
#   facet_grid(.~COG_category, scales = "free_x") +
#   
#   geom_bar(stat = "identity", position = position_dodge()) +
#   
#   theme_minimal() + 
#   
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

main_figure_COG_table <- df_plot_disease_asso_v2 %>%
  filter(COG_category == "M")

write.csv(main_figure_COG_table, "intermediate_files/main_figure_COG.csv", row.names = FALSE)

write.csv(df_plot_disease_asso_v2, "intermediate_files/supp_figure_COG.csv", row.names = FALSE)


