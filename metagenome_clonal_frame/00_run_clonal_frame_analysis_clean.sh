
# The configfile is a list of output directories, with each directory corresponding to a single putative sweep

CONFIGFILE=$1
OUTFOLDER_1=`awk '{if (NR=='${SLURM_ARRAY_TASK_ID}') print $1}' $CONFIGFILE`

#define directory where the script is called
HOMEFOLDER=`pwd`

#define and construct workspace on /tmp/$USER

THISWORKFOLDER=$TMPDIR
mkdir $THISWORKFOLDER/genomes
cp -R ${HOMEFOLDER}/clonal_frame_calling_clean $THISWORKFOLDER  

project_dir=$THISWORKFOLDER
OUTFOLDER=${HOMEFOLDER}/../Results/clonal_frame_analysis/VC_analysis/${OUTFOLDER_1}  #Edit me!
BASENAME=$OUTFOLDER_1  #Edit me!
export project_dir
export OUTFOLDER
export BASENAME

#make the output directory
mkdir -p $OUTFOLDER

# move to clonal_frame_calling
cd $THISWORKFOLDER/clonal_frame_calling_clean 
mkdir -p log

#Run prep code
sh 00_prepare_phybreak.sh
echo Done!

#Run first part of sweep pipeline
sh 01_run_phybreak.sh
echo Done!

#Run the tree building step
rm -R ${project_dir}/align/alignment_blocks
chgrp -R dome ${project_dir}
rsync -a --no-p ${project_dir}/./* $OUTFOLDER
rm -R $OUTFOLDER/align/alignment_blocks
module load R
cd $OUTFOLDER/clonal_frame_calling_clean   
Rscript fasta_to_phy.R $BASENAME
sbatch run_phyml_cf_gb.sh
