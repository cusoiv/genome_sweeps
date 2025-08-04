library(readr)
library(tidyr)
library(dplyr)

blast_names <- c("qseqid", "sseqid", "pident", "sstart", "send", 
                 "qstart", "qend", "evalue", "bitscore", "score", 
                 "qlen", "length")

SGB_vector <- read_csv("meta_data/SGB_list.csv") %>%
  pull(SGB)

blast_clonalframes_v_isolate_genomes_nucl_df <- read_table("raw_data/clonalframes_v_isolate_genomes_nucl.txt", col_names = FALSE)
colnames(blast_clonalframes_v_isolate_genomes_nucl_df) <- blast_names

proka_clonalframes_df <- read_delim("raw_data/clonalframes.tsv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

clonalframe_prokka_map_df <- proka_clonalframes_df %>%
  select(locus_tag, clonalframe) %>%
  rename(qseqid = locus_tag) %>%
  mutate(q_sweep_clonalframe = gsub("_cf_consensus", "", clonalframe))


########### Split up nucleotide blast results per SGB (maintain all subset genomes sequences)

for (SGB in SGB_vector){
  
  select_clonalframe_prokka_map_df <- clonalframe_prokka_map_df %>%
    filter(grepl(SGB, q_sweep_clonalframe))
  
  select_blast_clonalframes_v_isolate_genomes_nucl_df <- inner_join(blast_clonalframes_v_isolate_genomes_nucl_df, select_clonalframe_prokka_map_df)
  
  write.csv(select_blast_clonalframes_v_isolate_genomes_nucl_df, paste("subsetted_data/", SGB, "_clonalframes_v_isolate_genomes_nucl.csv", sep = ""), row.names = FALSE)
  
}

########### split up protein blast results per SGB

blast_clonalframes_all_v_all_prot_df <- read_table("raw_data/clonalframes_all_v_all_prot.txt", col_names = FALSE)
colnames(blast_clonalframes_all_v_all_prot_df) <- blast_names

clonalframe_prokka_map_q_df <- proka_clonalframes_df %>%
  select(locus_tag, clonalframe) %>%
  rename(qseqid = locus_tag) %>%
  mutate(q_sweep_clonalframe = gsub("_cf_consensus", "", clonalframe)) %>%
  select(-clonalframe)

clonalframe_prokka_map_s_df <- clonalframe_prokka_map_q_df %>%
  rename(sseqid = qseqid) %>%
  rename(s_sweep_clonalframe = q_sweep_clonalframe)

for (SGB in SGB_vector){
  
  select_clonalframe_prokka_map_q_df <- clonalframe_prokka_map_q_df %>%
    filter(grepl(SGB, q_sweep_clonalframe))
  
  select_clonalframe_prokka_map_s_df <- clonalframe_prokka_map_s_df %>%
    filter(grepl(SGB, s_sweep_clonalframe))
  
  select_blast_blast_clonalframes_all_v_all_prot_df <- inner_join(blast_clonalframes_all_v_all_prot_df, select_clonalframe_prokka_map_q_df) %>%
    inner_join(select_clonalframe_prokka_map_s_df)
  
  write.csv(select_blast_blast_clonalframes_all_v_all_prot_df, paste("subsetted_data/", SGB, "_clonalframes_all_v_all_prot.csv", sep = ""), row.names = FALSE)
  
}




