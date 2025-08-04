#define project_dir and file basename, need to be the same as that in phybreak_parameters.txt
output_dir=${OUTFOLDER}
project_dir=${project_dir}
genome_dir='/lisc/scratch/dome/xy43/UHGG_plus4/Data/VC_genomes/VC_O1_Shapiro_select'   #directory where your genomes are
input_contig_extension='fna.gz' #extension of genomes, can be zipped or unzipped
contig_extension='fna'   #contig_extension is the unzipped version of input_contig_extension
basename=${BASENAME}   #edit me! this is the basename of all your outputs
echo ${basename}.new

pop_infile_source='/lisc/scratch/dome/xy43/UHGG_plus4/Data/VC_genomes'  #edit me! this is where your genome input list is
pop_infile_name=${BASENAME}.txt #this is the file name of the genome input list

#use this if you want to select from genome quality file
#genome_info='/lisc/scratch/dome/xy43/UHGG_plus4/Results/drep_checkM_total/genomeInfo_total.csv' #this is where the N50 of the genomes are
#or you can define your own ref_iso, but it needs to be a complete genome or stitched together into a complete one
ref_iso='2010EL_1786'
ref_contig='2010EL_1786_1'



