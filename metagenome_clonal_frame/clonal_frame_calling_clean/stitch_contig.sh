configfile=./phybreak_config.sh
refisofile=./select_ref_isolate.sh
source ${configfile}
source ${refisofile}

echo $ref_iso
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < ${project_dir}/genomes/${ref_iso}.${contig_extension} > ${project_dir}/genomes/${ref_iso}_linear.fna
awk -F'\t' '{print $2}' ${project_dir}/genomes/${ref_iso}_linear.fna > ${project_dir}/genomes/newOH.fna
awk -v ref_contig="$ref_contig" 'BEGIN { ORS="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"; print ">" ref_contig "\n" } { print }' < ${project_dir}/genomes/newOH.fna > ${project_dir}/genomes/${ref_iso}_1.fna
rm ${project_dir}/genomes/${ref_iso}.fna
rm ${project_dir}/genomes/${ref_iso}_linear.fna
rm ${project_dir}/genomes/newOH.fna
mv ${project_dir}/genomes/${ref_iso}_1.fna ${project_dir}/genomes/${ref_iso}.${contig_extension}

