00_run_clonal_frame_analysis_clean.sh will call the contents in the folder clonal_frame_calling_clean to generate a consensus clonal frame for a putative sweep

01_run_bowtie_instrain_cf_MP1_MG.sh runs instrain and metaphlan for each metagenome against every clonal frame to calculate a a) metagenome to consensus clonal frame distance and b) determine if the metagenome contains only one strain for each SGB of interest

02_db_analysis_r2.R applies the 5X rule towards the results to find metagenome-isolate based sweeps

03_run_instrain_cf_compare.sh confirms all the metagenome-isolate sweeps with pairwise distances 
