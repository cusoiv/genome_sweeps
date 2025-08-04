#!/bin/bash

#SBATCH --job-name=coresimul
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --nice=0
#SBATCH --partition=basic
#SBATCH --output=log/%x-%A_%a.out
#SBATCH --error=log/%x-%A_%a.err
#SBATCH --time=06:00:00

# The configfile are the parameters to use
CONFIGFILE=$1
i=$((SLURM_ARRAY_TASK_ID+1))
rho=`awk -F',' '{if (NR=='${i}') print $1}' $CONFIGFILE`
delta=`awk -F',' '{if (NR=='${i}') print $2}' $CONFIGFILE`
scale=`awk -F',' '{if (NR=='${i}') print $3}' $CONFIGFILE`
coeff=`awk -F',' '{if (NR=='${i}') print $3}' $CONFIGFILE`

#define directory where the script is called
HOMEFOLDER=`pwd`

#define and construct workspace on /tmp/$USER
#mkdir /tmp/tmp_cluster
THISWORKFOLDER=$TMPDIR
cd $THISWORKFOLDER

#copy relevant files to workfolder 
cp ${HOMEFOLDER}/CoreSimul/*.py ${THISWORKFOLDER}
cp ${HOMEFOLDER}/CoreSimul/mycontrol_64_pop3.txt ${THISWORKFOLDER}
cp ${HOMEFOLDER}/CoreSimul/mytree_64_pop3.nwk ${THISWORKFOLDER}
cp ${HOMEFOLDER}/CoreSimul/ref.fa ${THISWORKFOLDER}

#edit control.txt according to input paramters
sed -i "s|.SLURM_ARRAY_TASK_ID|$SLURM_ARRAY_TASK_ID|g" mycontrol_64_pop3.txt
sed -i "s|.rho|$rho|g" mycontrol_64_pop3.txt
sed -i "s|.scale|$scale|g" mycontrol_64_pop3.txt
sed -i "s|.delta|$delta|g" mycontrol_64_pop3.txt
sed -i "s|.coeff|$coeff|g" mycontrol_64_pop3.txt

#run analysis
module load python3
python3 coresimul_master.py mycontrol_64_pop3.txt

#copy back to Results
mkdir -p ${HOMEFOLDER}/../Results/simulated_evolution/coresimul_64_pop3

rsync -aR ././test*/* ${HOMEFOLDER}/../Results/simulated_evolution/coresimul_64_pop3


