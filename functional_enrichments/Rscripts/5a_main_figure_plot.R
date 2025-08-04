library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)


main_figure_COG <- read_csv("intermediate_files/main_figure_COG.csv") %>%
  rename(Annotation_category = COG_category)

main_figure_PFAM <- read_csv("intermediate_files/main_figure_PFAM.csv") %>%
  rename(Annotation_category = PFAM_clean)

main_figure_df <- bind_rows(main_figure_COG, main_figure_PFAM) %>%
  mutate(Group = gsub("All_Sweeps", "All_Genes | All_Sweeps", Combined_ID)) %>%
  mutate(Group = gsub("Sweep_specific_Depleted", "Sweep_specific_genes | Depleted_Association", Group)) %>%
  mutate(Group = gsub("Sweep_specific_Enriched", "Sweep_specific_genes | Enriched_Association", Group)) %>%
  mutate(Group = gsub("Sweep_specific_Other", "Sweep_specific_genes | No_Association", Group)) %>%
  filter(Species != "s__Bacteroides_uniformis") %>%
  filter(Species != "s__Phocaeicola_vulgatus") %>%
  separate(Group, into = c("Gene_Group", "Sweep_Group"), sep = " \\| ", remove = FALSE)

main_figure_df$Annotation_category <- factor(main_figure_df$Annotation_category, levels = c("M", "Glyco_trans"))  

pdf(file = "main_figure_plot.pdf",   
    width = 8, 
    height = 4)

ggplot(data = main_figure_df, aes(x = Group, y = Percent, fill = Sweep_Group)) +
  
  facet_grid(Annotation_category~SGB*Gene_Group, space = "free_x", scales = "free_x") +
  
  geom_bar(stat = "identity", position = position_dodge()) +
  
  theme_minimal() + 
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(legend.position="top")

dev.off()