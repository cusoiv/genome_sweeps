#!/bin/bash
#
#SBATCH --job-name=phyml
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --nice=0
#SBATCH --partition=basic
#SBATCH --output=log/%x-%A_%a.out
#SBATCH --error=log/%x-%A_%a.err
#SBATCH --time=6-12:00:00

#usage
#sbatch -a 1-number of chunks ../Scripts/run_uclust_after_prodigal_rep_seq.sh <CONFIGFILE>

#The configfile is a list of folder files
#CONFIGFILE=$1
#FOLDER=`awk '{if (NR=='${SLURM_ARRAY_TASK_ID}') print $1}' $CONFIGFILE`
#BASENAME=`echo $FOLDER | awk  'BEGIN{FS="[.-]"; OFS=""}; {print $1,$2,"_",$NF}'`
#echo $BASENAME

##define directory where the script is called
HOMEFOLDER=`pwd`
#BASENAME=clonal_anicluster90diffH  #edit me!

#define and construct workspace on /tmp/$USER
#rm -R /tmp/tmp_2
#mkdir /tmp/tmp_2
#THISWORKFOLDER=/tmp/tmp_2
THISWORKFOLDER=$TMPDIR
#mkdir $THISWORKFOLDER/output

#copy core fasta to $TMPDIR
cp ${OUTFOLDER}/align/${BASENAME}.core.fasta $THISWORKFOLDER
cp ${OUTFOLDER}/align/${BASENAME}.core.fasta.phy $THISWORKFOLDER
cp ${OUTFOLDER}/align/* $THISWORKFOLDER
cp -R ${OUTFOLDER}/clonal_frame_calling_clean $THISWORKFOLDER

#run phyml
module load phyml
cd $THISWORKFOLDER
phyml -i ${BASENAME}.core.fasta.phy -m 'GTR' -t 'e' -a 'e' -f 'm' -v 'e' > ${BASENAME}.core.phy_phyml_stat.txt #GTR+G+I, see https://lorenzogatti.me/2016_FiPS_Tutorials/solutions_tutorial01.html

#calculate tv/ts rate 
AC=`grep 'A <-> C' ${BASENAME}.core.fasta.phy_phyml_stats.txt | awk '{print $NF}'`
AG=`grep 'A <-> G' ${BASENAME}.core.fasta.phy_phyml_stats.txt | awk '{print $NF}'`
AT=`grep 'A <-> T' ${BASENAME}.core.fasta.phy_phyml_stats.txt | awk '{print $NF}'`
CG=`grep 'C <-> G' ${BASENAME}.core.fasta.phy_phyml_stats.txt | awk '{print $NF}'`
CT=`grep 'C <-> T' ${BASENAME}.core.fasta.phy_phyml_stats.txt | awk '{print $NF}'`
GT=`grep 'G <-> T' ${BASENAME}.core.fasta.phy_phyml_stats.txt | awk '{print $NF}'`

kappa=`awk "BEGIN {z=2*($AG+$CT)/($AC+$AT+$CG+$GT); print z}"`


#run ClonalframeML
module load conda
conda activate clonalframeml
ClonalFrameML ${BASENAME}.core.fasta.phy_phyml_tree.txt ${BASENAME}.core.fasta ${BASENAME}.output -kappa ${kappa} -emsim 100  > ${BASENAME}.log.txt

#run Gubbins
#module load gubbins
#module load raxml
#run_gubbins.py ${BASENAME}.core.fasta -s ${BASENAME}.core.fasta.phy_phyml_tree.txt

#produce masked alignments after clonalframeml & gubbins
conda activate /lisc/user/yu/envs/maskrc-svg
~/maskrc-svg.py --aln ${BASENAME}.core.fasta --out ${BASENAME}.core.cf.maskrc.aln ${BASENAME}.output   #clonalframeml
#maskrc-svg.py --aln ${BASENAME}.core.fasta --out ${BASENAME}.core.gb.maskrc.aln --gubbins ${BASENAME}.core  #gubbins

#convert fasta to phy format
#clonalframeml
awk -v RS=">" '{gsub(/\n/," ");} NF' ${BASENAME}.core.cf.maskrc.aln > test
cat <( head -1 ${BASENAME}.core.fasta.phy ) test > ${BASENAME}.core.cf.maskrc.phy

#gubbins
#awk -v RS=">" '{gsub(/\n/," ");} NF' ${BASENAME}.core.gb.maskrc.aln > test
#cat <( head -1 ${BASENAME}.core.fasta.phy ) test > ${BASENAME}.core.gb.maskrc.phy

#build phyml tree from clonalframe and gubbins results
phyml -i ${BASENAME}.core.cf.maskrc.phy -m 'GTR' -t 'e' -a 'e' -f 'm' -v 'e' > ${BASENAME}.core.cf.maskrc.phy_phyml_stat.txt
#phyml -i ${BASENAME}.core.gb.maskrc.phy -m 'GTR' -t 'e' -a 'e' -f 'm' -v 'e' > ${BASENAME}.core.gb.maskrc.phy_phyml_stat.txt

#get consensus sequence for clonalframe
module load python3
cd $THISWORKFOLDER/clonal_frame_calling_clean
sed "s|.project_dir|$THISWORKFOLDER|g" phybreak_parameters_template.txt > phybreak_parameters.txt
sed -i "s|.basename|$BASENAME|g" phybreak_parameters.txt
python3 phybreak8.core_to_cf.py

#rsync back results
chgrp -R $THISWORKFOLDER
rsync -a --no-p $THISWORKFOLDER/* $OUTFOLDER/align/
 
