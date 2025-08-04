library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)


supp_figure_COG <- read_csv("intermediate_files/supp_figure_COG.csv") %>%
  rename(Annotation_category = COG_category)

supp_figure_PFAM <- read_csv("intermediate_files/supp_figure_PFAM.csv") %>%
  rename(Annotation_category = PFAM_clean)

supp_figure_df <- bind_rows(supp_figure_COG, supp_figure_PFAM) %>%
  mutate(Group = gsub("All_Sweeps", "All_Genes | All_Sweeps", Combined_ID)) %>%
  mutate(Group = gsub("Sweep_specific_Depleted", "Sweep_specific_genes | Depleted_Association", Group)) %>%
  mutate(Group = gsub("Sweep_specific_Enriched", "Sweep_specific_genes | Enriched_Association", Group)) %>%
  mutate(Group = gsub("Sweep_specific_Other", "Sweep_specific_genes | No_Association", Group)) %>%
  separate(Group, into = c("Gene_Group", "Sweep_Group"), sep = " \\| ", remove = FALSE)

supp_figure_df$Annotation_category <- factor(supp_figure_df$Annotation_category, levels = c("M", "L", "Glyco_trans", "Phage_integrase", "DDE_Tnp", "Methylase_S"))  

pdf(file = "supp_figure_plot.pdf",   
    width = 12, 
    height = 16)

ggplot(data = supp_figure_df, aes(x = Group, y = Percent, fill = Sweep_Group)) +
  
  facet_grid(Annotation_category~SGB*Gene_Group, space = "free_x", scales = "free_x") +
  
  geom_bar(stat = "identity", position = position_dodge()) +
  
  theme_minimal() + 
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(legend.position="top")

dev.off()
