#!/bin/bash
#usage
#sbatch -a 1-number of chunks run_bowtie_instrain_cf_MP1_MG.sh $CONFIGFILE #CONFIGFILE are lists of metagenomes.txt, with the file containing NCBI accession numbers corresponding to the metagenome
#define directory where the script is called
HOMEFOLDER=`pwd`
CONFIGFILE=$1

#define and construct workspace on /tmp/$USER
THISWORKFOLDER=$TMPDIR
#mkdir /tmp/tmp_${SLURM_ARRAY_TASK_ID}
#THISWORKFOLDER=/tmp/tmp_${SLURM_ARRAY_TASK_ID}
GENOMENAME=`awk '{if (NR=='${SLURM_ARRAY_TASK_ID}') print $1}' $CONFIGFILE`
cp $HOMEFOLDER/../Results/clonal_groups/metagenome_lists/${GENOMENAME}.txt $THISWORKFOLDER
mn=`grep -n $GENOMENAME $HOMEFOLDER/../Results/clonal_groups/metagenome_lists/MG.list | cut -d':' -f1`
echo $mn
module load sratoolkit/
cd $THISWORKFOLDER

#download metagenome of interest

while read line; do
	prefetch $line
	fasterq-dump $line
done<${GENOMENAME}.txt

cat *_1.fastq > ${GENOMENAME}_1.fastq
cat *_2.fastq > ${GENOMENAME}_2.fastq

#Run instrain --somehow metaphlan modifies original fastq so need to run instrain first!
module load bowtie2  #version 2.5.1
module load conda
module load samtools
conda activate instrain-1.7.5

for i in {1..6};do
#for i in 6; do
	basename=consensus_db_${i}
	rsync -a $HOMEFOLDER/../Results/clonal_frame_analysis/clonal_SGBdiffH_sweep/new_consensus_cfs/db_${i}/${basename}* $THISWORKFOLDER/
	bowtie2 -X 2000 -x ${basename} -1 ${GENOMENAME}_1.fastq -2 ${GENOMENAME}_2.fastq -S ${GENOMENAME}.sam --no-mixed --very-sensitive -p 8
	inStrain profile ${GENOMENAME}.sam ${basename}.fna -g ${basename}.genes.fna -o ${GENOMENAME}.IS -p 8 -s ${basename}.stb
	cd ${GENOMENAME}.IS
	tar czf raw_data.tar.gz raw_data --remove-files
	cd ..
#	fn=$(echo "scale=0; $SLURM_ARRAY_TASK_ID / 100 + 1" | bc -l)
	fn=$((mn / 100 + 1))
	echo $fn
	mkdir -p $HOMEFOLDER/../Results/instrain/db_${i}_MQ1/MG${fn}
	rsync -azR  ./*.IS  $HOMEFOLDER/../Results/instrain/db_${i}_MQ1/MG${fn}
done

#Run metaphlan
conda activate metaphlan-4.0.6
mkdir -p $THISWORKFOLDER/marker_dir
rsync -aR $HOMEFOLDER/../../strainphlan_test/Data/meta4_beta1_vJan21/./* $THISWORKFOLDER/marker_dir
suffix=marker
metaphlan ${GENOMENAME}_1.fastq ${GENOMENAME}_2.fastq --input_type fastq -s ${GENOMENAME}_${suffix}.sam.bz2 --bowtie2out ${GENOMENAME}_${suffix}.bowtie2.bz2 -o ${GENOMENAME}_${suffix}_profile.tsv --bowtie2db marker_dir --index mpa_vJan21_CHOCOPhlAnSGB_202103 --nproc ${SLURM_CPUS_PER_TASK} --min_mapq_val 20
metaphlan ${GENOMENAME}_${suffix}.bowtie2.bz2 --input_type bowtie2out -o ${GENOMENAME}_${suffix}_profile.tsv --bowtie2db marker_dir --index mpa_vJan21_CHOCOPhlAnSGB_202103 --nproc ${SLURM_CPUS_PER_TASK} --min_mapq_val 20 --tax_lev 't'

#sync results to new folder
mkdir -p $HOMEFOLDER/../Results/strainphlan/MG${fn}/${GENOMENAME}
cp *sam.bz2 $HOMEFOLDER/../Results/strainphlan/MG${fn}/${GENOMENAME}
cp *_profile.tsv $HOMEFOLDER/../Results/strainphlan/MG${fn}/${GENOMENAME}

