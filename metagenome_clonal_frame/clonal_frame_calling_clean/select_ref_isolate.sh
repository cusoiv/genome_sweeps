configfile=./phybreak_config.sh
source ${configfile}
#genome_info='/scratch/xy43/UHGG_BN10/Results/drep_checkM_total/genomeInfo_total.csv'

cut -d$'\t' -f1 ${project_dir}/${pop_infile_name} > pop_names

#find ref_isolate
grep -f pop_names $genome_info > genome_info_select

ref_iso=`basename -s .${contig_extension} $(awk -F',' -v max=0 '{if($NF>max){want=$1; max=$NF}}END{print want}' genome_info_select)`
ref_contig=${ref_iso}_1

rm pop_names
rm genome_info_select


 
