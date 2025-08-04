module load conda
conda activate /lisc/scratch/dome/xy43/miniconda/envs/popcogent_test
configfile=./phybreak_config.sh
source ${configfile}
module load perl
python phybreak1.generate_maf.py
python phybreak2.maf_to_fasta.py

rsync -aR --no-p ${project_dir}/./* ${output_dir}



