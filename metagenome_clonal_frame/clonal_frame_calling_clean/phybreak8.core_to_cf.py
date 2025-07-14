# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 14:56:07 2021

@author: Xiaoqian
"""


import os
#os.chdir("/lisc/scratch/dome/xy43/UHGG_plus4/Results/clonal_frame_analysis/VC_analysis/VC_haiti/clonal_frame_calling_clean")

## Collect parameters
project_dir = ""
input_contig_dir = ""
contig_dir = ""
contig_extension = ""
output_prefix = ""
pop_infile_name = ""
ref_iso = ""
ref_contig = ""
focus_population = ""
len_block_threshold = 0
gap_prop_thresh = 0.0
window_size = 0
overlap = 0
MUGSY_source = ""
phyML_loc = ""
phyML_properties = ""
ape_loc = ""
percentile_threshold = 0.0
min_physplit_window_size = 0


parameter_file = open("phybreak_parameters.txt","r")
for line in parameter_file:
	line = line.strip().split(" = ")
	if len(line) > 1:
		if line[0] == "project_dir":
			project_dir = line[1].split(" #")[0]
		elif line[0] == "input_contig_dir":
			input_contig_dir = line[1].split(" #")[0]
		elif line[0] == "contig_dir":
			contig_dir = line[1].split(" #")[0]
		elif line[0] == "input_contig_extension":
			contig_extension = line[1].split(" #")[0]
		elif line[0] == "output_prefix":
			output_prefix = line[1].split(" #")[0]
		elif line[0] == "pop_infile_name":
			pop_infile_name = line[1].split(" #")[0]
		elif line[0] == "ref_iso":
			ref_iso = line[1].split(" #")[0]
		elif line[0] == "ref_contig":
			ref_contig = line[1].split(" #")[0]
		elif line[0] == "focus_population":
			focus_population = line[1].split(" #")[0]
		elif line[0] == "len_block_threshold":
			len_block_threshold = int(line[1].split(" #")[0])
		elif line[0] == "gap_prop_thresh":
			gap_prop_thresh = float(line[1].split(" #")[0])
		elif line[0] == "window_size":
			window_size = int(line[1].split(" #")[0])
		elif line[0] == "window_overlap":
			overlap = int(line[1].split(" #")[0])
		elif line[0] == "MUGSY_source":
			MUGSY_source = line[1].split(" #")[0]
		elif line[0] == "phyML_loc":
			phyML_loc = line[1].split(" #")[0]
		elif line[0] == "phyML_properties":
			phyML_properties = line[1].split(" #")[0]
		elif line[0] == "ape_loc":
			ape_loc = line[1].split(" #")[0]
		elif line[0] == "percentile_threshold":
			percentile_threshold = float(line[1].split(" #")[0])
		elif line[0] == "min_physplit_window_size":
			min_physplit_window_size = int(line[1].split(" #")[0])
parameter_file.close()

#these directories will generate if they do not already exist
input_dir = project_dir+"align/"
alignment_dir = input_dir+"alignment_blocks/"
phy_split_dir = input_dir+"phy_split/"
tree_dir = input_dir+"trees/"
msa_out_dir = input_dir+"phybreak_blocks/"

#these are output file names
strain_list_filename = "strain_names.txt"
MSA_name = output_prefix+".core.fasta"
MSA_name_cf_maskrc = output_prefix+".core.cf.maskrc.aln"
consensus_core = output_prefix+".core.cf.consensus.fasta"



def find_consensus_seq(msa_dict):
	nt_dict = {}
	

	for strain in msa_dict:
#		print (strain)
		seq = msa_dict[strain]
		for i in range(0,len(seq)):

			try:
				nt_dict[i]  #if nt_dict[i] exists
			except:
				nt_dict[i] = {'N':0,'?':0,'A':0,'T':0,'G':0,'C':0}
			nt = seq[i]
			if nt == '?':
				nt_dict[i]['?'] += 1
			elif nt == 'N':
				nt_dict[i]['N'] += 1
			elif nt == 'A':
				nt_dict[i]['A'] += 1
			elif nt == 'T':
				nt_dict[i]['T'] += 1
			elif nt == 'G':
				nt_dict[i]['G'] += 1
			elif nt == 'C':
				nt_dict[i]['C'] += 1
	consensus_seq=""		
	for i in nt_dict:
		sn=0
		for l in 'ATGC':
			sn=sn+nt_dict[i][l]
		if sn==sum(nt_dict[i].values()):#When there is no recombination in any of the genomes involved
			consensus_seq=consensus_seq+max(nt_dict[i],key=nt_dict[i].get)
	
	return consensus_seq



				
msa = open(input_dir+MSA_name_cf_maskrc,"r")
#msa=open('testmaf.core.keepgap.fasta')
seq_dict = {}
head = ""
for line in msa:
	line = line.strip()
	if line[0] == ">":
		head = line[1:len(line)]
	else:
		try:
			seq_dict[head] += line
		except:
			seq_dict[head] = line
msa.close()		

consensus_seq=find_consensus_seq(seq_dict)

outfile=open(input_dir+consensus_core,"w")
#outfile=open("test_consensus_core.fasta","w")
outfile.write(">consensus_core_cf\n"+consensus_seq+"\n")
outfile.close()
