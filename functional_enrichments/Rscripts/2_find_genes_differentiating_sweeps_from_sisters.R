library(readr)
library(tidyr)
library(dplyr)

sweep_genome_relationships <- read_csv("meta_data/cleaned_sweep_metadata.csv") %>%
  rename(s_genome = genome) %>%
  rename(s_sweep_relationship = sweep) %>%
  rename(relationship_type = relationship) %>%
  
  group_by(s_sweep_relationship, relationship_type) %>%
  mutate(n_genomes = length(unique(s_genome))) %>%
  ungroup()

prokka_isolate_genomes_df <- read_delim("raw_data/isolate_genomes.tsv", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)

isolate_genomes_prokka_map_df <- prokka_isolate_genomes_df  %>%
  select(locus_tag, genome) %>%
  rename(sseqid = locus_tag) %>%
  rename(s_genome = genome)

emapper_annotations_cf_df <- read_delim("raw_data/clonalframes.emapper.annotations", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE, skip = 4) 
colnames(emapper_annotations_cf_df)[1] <- "qseqid"

######

input_dir <- "subsetted_data/"

files <- list.files(input_dir, pattern = "_nucl.csv")

df_sweep_specific_list <- list()

i <- 1

for (file in files){
  
    SGB_blast <- read_csv(paste(input_dir, file, sep = ""))
    
    SGB_compiled_comp_df <- inner_join(SGB_blast, isolate_genomes_prokka_map_df) %>%
      inner_join(sweep_genome_relationships, relationship = "many-to-many") %>%
      
      # very import filter so that we are only comparing blast results within sweep and to sisters of sweep
      filter(q_sweep_clonalframe == s_sweep_relationship) %>%
      
      # clean up
      mutate(paln = (length / qlen)*100) %>%
      select(qseqid, sseqid, paln, pident, q_sweep_clonalframe, s_genome, n_genomes, s_sweep_relationship, relationship_type)
    
    # select blast hits with identical sequences between clonalframe and all
    
    gene_conserved_within_sweep <- SGB_compiled_comp_df %>%
      
      filter(relationship_type == "within_sweep") %>%
      
      filter(paln == 100) %>%
      filter(pident == 100) %>%
      
      group_by(qseqid, q_sweep_clonalframe, s_sweep_relationship) %>%
      mutate(n_hit_genomes = length(unique(s_genome))) %>%
      ungroup() %>%
      
      filter(n_hit_genomes == n_genomes) %>%
      select(qseqid) %>%
      distinct() %>%
      pull(qseqid)
    
    
    # select for genes that are missing in the siter genomes (look at some different cutoffs)
    # first maybe just look at some distrubions of hits to sister genomes
    
    genes_found_in_closest_sisters <- SGB_compiled_comp_df %>%
      
      filter(relationship_type == "closest_sister") %>%
      
      # if alignment identity below 30%, gene is considered not found
      
      filter(paln > 30) %>%
      
      pull(qseqid)
    
    
    # intersection of genes conserved with sweep, but not found in sister
    sweep_specific_genes <- gene_conserved_within_sweep[!(gene_conserved_within_sweep %in% genes_found_in_closest_sisters)]
    
    # annotations of sweep specific genes
    sweep_specific_genes_annotations <- emapper_annotations_cf_df %>%
      filter(qseqid %in% sweep_specific_genes)
    
    df_sweep_specific_list[[i]] <- sweep_specific_genes_annotations
    i <- i + 1
    
}

map_SGB_df <- read_delim("raw_data/clonalframes.tsv", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE) %>%
  
  select(locus_tag, gene, product, clonalframe) %>%
  
  separate(clonalframe, into = c("SGB", "sweep_type", "sweep_num"), remove = FALSE) %>%
  
  mutate(Sweep = gsub("_cf_consensus", "", clonalframe)) %>%
  
  rename(qseqid = locus_tag) %>%
  rename(prokka_gene = gene) %>%
  rename(prokka_annotation = product)
  
df_sweep_specific_compiled <- bind_rows(df_sweep_specific_list) %>%
  inner_join(map_SGB_df) %>%
  select(-seed_ortholog, -evalue, -score, -eggNOG_OGs, -GOs, -EC) 

df_sweep_specific_compiled$SGB <- gsub("diffH", "", df_sweep_specific_compiled$SGB)
df_sweep_specific_compiled$Sweep <- gsub("diffH", "", df_sweep_specific_compiled$Sweep)


write.csv(df_sweep_specific_compiled, "intermediate_files/genes_differentiating_sweeps_from_sisters.csv", row.names = FALSE)

