library(readr)
library(tidyr)
library(dplyr)

input_dir <- "subsetted_data/"

prot_blast_list <- list.files(input_dir, pattern = "_prot.csv")

df_prot_sweep_specific_list <- list()

i <- 1

for (file in prot_blast_list){

  blast_cf_v_cf_prot_df <- read_csv(paste(input_dir, file, sep = ""))
  
  n_sweeps <- length(unique(c(blast_cf_v_cf_prot_df$q_sweep_clonalframe, blast_cf_v_cf_prot_df$s_sweep_clonalframe)))
  
  if(n_sweeps >= 3){
  
    cf_specific_prots_df <- blast_cf_v_cf_prot_df %>%
      
      mutate(paln = (length / qlen)*100) %>%
      
      select(qseqid, sseqid, paln, pident, q_sweep_clonalframe, s_sweep_clonalframe) %>%
      
      filter(paln > 60) %>%
      filter(pident > 60) %>%
   
      group_by(qseqid) %>%
      mutate(n_sweep_hits = length(unique(s_sweep_clonalframe))) %>%
      ungroup() %>%
      
      filter(n_sweep_hits == 1) %>%
      filter(q_sweep_clonalframe == s_sweep_clonalframe) %>%
      
      filter(pident == 100) %>%
      filter(paln == 100) %>%
      
      select(qseqid, q_sweep_clonalframe) %>%
      
      separate(qseqid, into=c("rm", "ORF_number"), sep = "_", remove = FALSE) %>%
      
      select(-rm) %>%
      
      mutate(ORF_number = as.numeric(ORF_number)) 
    
      
      df_prot_sweep_specific_list[[i]] <- cf_specific_prots_df
      i <- i + 1
  }
}

compiled_prot_sweep_specific_df <- bind_rows(df_prot_sweep_specific_list)


genome_sweep_map <- read_csv("meta_data/cleaned_sweep_metadata.csv") 

genome_sweep_map <- genome_sweep_map %>%
  
  group_by(sweep) %>%
  mutate(n_genomes = length(unique(genome))) %>%
  ungroup() %>%
  
  select(-relationship)

prokka_annotations_cf_df <- read_delim("raw_data/clonalframes.tsv", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)

prokka_annotations_genomes_df <- read_delim("raw_data/isolate_genomes.tsv", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE) %>%
  
  inner_join(genome_sweep_map)

emapper_annotations_cf_df <- read_delim("raw_data/clonalframes.emapper.annotations", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE, skip = 4) 
colnames(emapper_annotations_cf_df)[1] <- "locus_tag"


final_results_clean <- compiled_prot_sweep_specific_df %>%
  
  rename(locus_tag = qseqid) %>%
  
  inner_join(prokka_annotations_cf_df) %>%
  
  left_join(emapper_annotations_cf_df) %>%
  
  # clean up
  
  select(-seed_ortholog, -evalue, -score, -eggNOG_OGs, -max_annot_lvl, -GOs, -EC, -KEGG_ko, -KEGG_Pathway, -KEGG_Reaction, -KEGG_rclass, -BRITE, -KEGG_TC, -BiGG_Reaction) %>%
  select(locus_tag, clonalframe, ORF_number, length_bp, COG_category, Preferred_name, gene,  Description, product,  PFAMs) %>%
  
  rename(Clonalframe = clonalframe) %>%
  rename(ORF_id = locus_tag) %>%
  rename(Gene_length = length_bp) %>%
  rename(Gene_name_EggNOG = Preferred_name) %>%
  rename(Gene_name_Prokka = gene) %>%
  rename(Description_EggNOG = Description) %>%
  rename(Description_Prokka = product) %>% 
  
  mutate_all(~ifelse(is.na(.), "-", .))


write.csv(final_results_clean, "intermediate_files/unique_prots_across_sweep.csv", row.names = FALSE)


