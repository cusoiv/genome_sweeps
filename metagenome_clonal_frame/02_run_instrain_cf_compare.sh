#!/bin/bash
#usage
#sbatch -a 1-number of chunks run_instrain_cf_compare.sh $CONFIGFILE #CONFIGFILE is a list of putative sweeps
#define directory where the script is called
HOMEFOLDER=`pwd`
CONFIGFILE=$1

#define and construct workspace on /tmp/$USER
THISWORKFOLDER=$TMPDIR
#mkdir /tmp/tmp_${SLURM_ARRAY_TASK_ID}
#THISWORKFOLDER=/tmp/tmp_${SLURM_ARRAY_TASK_ID}
GENOMENAME_1=`awk '{if (NR=='${SLURM_ARRAY_TASK_ID}') print $1}' $CONFIGFILE`
GENOMENAME=`echo $GENOMENAME_1 | sed 's/_consensus.txt//g'`
CF=`echo $GENOMENAME | sed 's|diffH_|diffH/|g'`

echo $GENOMENAME_1
echo $GENOMENAME
echo $CF

#copy relevant IS folders
#IS list
cp $HOMEFOLDER/../Results/clonal_groups/compare_select_lists/${GENOMENAME_1} $THISWORKFOLDER
cp $HOMEFOLDER/generate_pairs.sh $THISWORKFOLDER
#move to working directory
cd $THISWORKFOLDER

#find which db the consensus belongs to 
cp $HOMEFOLDER/../Results/clonal_frame_analysis/clonal_SGBdiffH_sweep/new_consensus_cfs/SGB_cf_consensus.list .
gf=`grep ${GENOMENAME}_ SGB_cf_consensus.list`
db=`dirname $gf`

ls -d $HOMEFOLDER/../Results/instrain/${db}_MQ1/MG*/./*.IS > list1
ls -d $HOMEFOLDER/../Results/instrain/${db}_MQ1/${CF}/./*.IS > list2
ls -d $HOMEFOLDER/../Results/instrain/${db}_MQ1/${CF}_sister/./*.IS > list3

cat list1 list2 list3 > ISlist_all
grep -f ${GENOMENAME_1} ISlist_all > IS_select.list
while read line;do
	rsync -aR $line $THISWORKFOLDER 
done<IS_select.list

for fIS in *.IS; do
        cd $fIS
        tar -xf raw_data.tar.gz
        cd ..
done

module load conda
conda activate instrain-1.7.5

basename=consensus_${db}
rsync -a $HOMEFOLDER/../Results/clonal_frame_analysis/clonal_SGBdiffH_sweep/new_consensus_cfs/${db}/${basename}* $THISWORKFOLDER
sh generate_pairs.sh

while read line; do
        cn=`echo $line | sed 's/ /_/g'`
	echo $cn
        inStrain compare -i $line -o ${cn}_compare -s ${basename}.stb -p 1 --genome ${GENOMENAME}_cf_consensus.fasta --skip_plot_generation
done<pairs.list


ls *_compare/output/*_compare_genomeWide_compare.tsv> compare.list

f1=`head -1 compare.list`
head -1 $f1 >> instrainComparer_genomeWide_compare.tsv
while read line; do
        tail -n+2 $line >> instrainComparer_genomeWide_compare.tsv
done<compare.list

#sync results to new folder
mkdir -p $HOMEFOLDER/../Results/instrain/${db}_MQ1/${CF}_compare
cp instrainComparer_genomeWide_compare.tsv $HOMEFOLDER/../Results/instrain/${db}_MQ1/${CF}_compare


