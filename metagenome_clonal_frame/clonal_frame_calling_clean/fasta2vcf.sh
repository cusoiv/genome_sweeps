module load conda
conda activate /home/user/yu/my-conda-envs/snp-sites
for f in *.fasta; do
	fn=`basename -s .fasta $f`
	snp-sites -v -o ${fn}.vcf ${fn}.fasta
done
