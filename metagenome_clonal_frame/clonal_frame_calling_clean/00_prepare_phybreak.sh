configfile=./phybreak_config.sh
source ${configfile}


#copy genomes from population to project_dir/genomes
echo ${pop_infile_source}/${pop_infile_name}
mkdir -p ${project_dir}
cp ${pop_infile_source}/${pop_infile_name} ${project_dir}/
mkdir -p ${project_dir}/genomes
cut -d$'\t' -f1 ${project_dir}/${pop_infile_name} > pop_names

#pop_names is now a list of genomes
while read line;do
 echo $line
 cp ${genome_dir}/${line}.${input_contig_extension} ${project_dir}/genomes
done<pop_names

#unzip
gunzip ${project_dir}/genomes/*.${input_contig_extension}
rm pop_names

#stitch ref_genome together
#bash stitch_contig.sh

#import ref_iso if you need to select from genome quality file
#refisofile=./select_ref_isolate.sh
#source ${refisofile}
#you can also define your own ref_iso if you want without the selection process, see config file

#edit phybreak_parameter file

sed "s|.project_dir|$project_dir|g" phybreak_parameters_template.txt > phybreak_parameters.txt
sed -i "s|.pop_infile_name|$pop_infile_name|g" phybreak_parameters.txt
sed -i "s|.basename|$basename|g" phybreak_parameters.txt
sed -i "s|.focus_population|$focus_population|g" phybreak_parameters.txt
sed -i "s|.ref_iso|$ref_iso|g" phybreak_parameters.txt
sed -i "s|.ref_contig|$ref_contig|g" phybreak_parameters.txt
sed -i "s|.total_jobs|$total_jobs|g" phybreak_parameters.txt


