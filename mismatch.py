#!/usr/bin/env python
# Usage: python mismatch.py -q SAM_input
# This script is used for Jasleen et al. (2019) paper
# More than 1 SAM files can be put in the argument, like: python mismatch.py -q SAM_file1 SAM_file2
# Alignment file (SAM format) for both single-end reads and paired-end reads can be used
# This script count the number of mismatches added to the 3' end of reads
# Output will contain: Counts for 3'-end mismatches; 
#                      nucleotide pereferences for 3'-end nucleotides; 
#                      percentage of mismatch at every position;
#                      strandedness percentage of mismatch at every position  
# Please contact Feng Wang (wangfeng@iu.edu) if you have any questions

from __future__ import division
from Bio.Seq import Seq
import re
import argparse

def terminal_misincorporation(file):
	with open(file,'rb') as f:
		forward_3t_counter = {0:0,1:0,2:0,3:0}
		reverse_3t_counter = {0:0,1:0,2:0,3:0}
		forward_3t_NT = {}
		reverse_3t_NT = {}
		mismatch_pos = {"forward":{}, "reverse":{}}
		total_pos = {"forward":{}, "reverse":{}}
		
		for line in f:
			if not line.startswith("@"):
				counter = 0
				alignment = line.strip().split('\t')
				ID = alignment[0]
				FLAG = alignment[1]
				Read_length = len(alignment[9])
				PE_size = int(alignment[8])
				Read = alignment[9]
				
				if FLAG == "0":	
					MDZ = alignment[12].split(":")[-1]
					Mismatch = map(int, re.split(r'[ACGTN]',MDZ))
					if Read_length not in mismatch_pos["forward"]:
						mismatch_pos["forward"][Read_length] = {}
						for i in range(1, Read_length + 1):
							mismatch_pos["forward"][Read_length][i] = 0
					if Read_length not in total_pos["forward"]:
						total_pos["forward"][Read_length] = 1
					else:
						total_pos["forward"][Read_length] += 1
					
					pos = 0
					for interval in Mismatch[:-1]:
						pos += interval + 1
						mismatch_pos["forward"][Read_length][pos] += 1
					while Mismatch.pop() == 0:
						counter += 1
					forward_3t_counter[counter] += 1
					sequence = Read
					MisNT = sequence[0-counter:]
					if counter >0:
						if MisNT not in forward_3t_NT:
							forward_3t_NT[MisNT] = 1
						else:
							forward_3t_NT[MisNT] += 1
				
				elif FLAG == "16":
					MDZ = alignment[12].split(":")[-1]
					Mismatch = map(int, re.split(r'[ACGTN]',MDZ))
					if Read_length not in mismatch_pos["reverse"]:
						mismatch_pos["reverse"][Read_length] = {}
						for i in range(1, Read_length + 1):
							mismatch_pos["reverse"][Read_length][i] = 0
					if Read_length not in total_pos["reverse"]:
						total_pos["reverse"][Read_length] = 1
					else:
						total_pos["reverse"][Read_length] += 1
					
					pos = 0
					for interval in Mismatch[::-1][:-1]:
						pos += interval + 1
						mismatch_pos["reverse"][Read_length][pos] += 1
					while Mismatch.pop(0) == 0:
						counter += 1
					reverse_3t_counter[counter] += 1
					sequence = Seq(Read).reverse_complement()
					MisNT = sequence[0-counter:]
					if counter > 0:
						if MisNT not in reverse_3t_NT:
							reverse_3t_NT[MisNT] = 1
						else:
							reverse_3t_NT[MisNT] += 1
				
				elif FLAG == "99":
					first_sequence = Read
					first_MDZ = alignment[12].split(":")[-1]
					first_Mismatch = map(int, re.split(r'[ACGTN]',first_MDZ))
					next_alignment = f.next().strip().split('\t')
					next_sequence = Seq(next_alignment[9])
					next_MDZ = next_alignment[12].split(':')[-1]
					next_Mismatch = map(int, re.split(r'[ACGTN]',next_MDZ))
					next_ID = next_alignment[0]
					next_FLAG = next_alignment[1]
					
					# count reads containing 3' mismatch
					for i in next_Mismatch[::-1]:
						if i ==0:
							counter += 1
						else:
							break
					forward_3t_counter[counter] += 1
					sequence = Read
					MisNT = sequence[0-counter:]
					if counter >0:
						if MisNT not in forward_3t_NT:
							forward_3t_NT[MisNT] = 1
						else:
							forward_3t_NT[MisNT] += 1

					# Begin counting mismatches at different positions
					if next_ID == ID and next_FLAG == "147":
						if PE_size <= 160:
							head_length = tail_length = PE_size - 80
							shared_length = 160 - PE_size
							if first_sequence[80 - shared_length:80] == next_sequence[0:shared_length]:
								if PE_size not in mismatch_pos["forward"]:
									mismatch_pos["forward"][PE_size] = {}
									for i in range(1, PE_size + 1):
										mismatch_pos["forward"][PE_size][i] = 0
								if PE_size not in total_pos["forward"]:
									total_pos["forward"][PE_size] = 1
								else:
									total_pos["forward"][PE_size] += 1
								
								next_Mismatch[0] += head_length
								pos = 0
								adj_next_Mismatch = []
								for interval in next_Mismatch:
									pos += interval + 1
									if pos > 80:
										adj_next_Mismatch.append(interval)
								pos = 0
								for interval in first_Mismatch:
									pre_pos = pos
									pos += interval + 1
									if pos > head_length:
										adj_next_Mismatch[0] -= pre_pos
										break
								Mismatch = first_Mismatch[:-1] + adj_next_Mismatch
								pos = 0
								for interval in Mismatch[:-1]:
									pos += interval + 1
									mismatch_pos["forward"][PE_size][pos] += 1
							else:
								continue
						else:
							gap = PE_size - 160
							if PE_size not in mismatch_pos["forward"]:
								mismatch_pos["forward"][PE_size] = {}
								for i in range(1, PE_size + 1):
									mismatch_pos["forward"][PE_size][i] = 0
							if PE_size not in total_pos["forward"]:
								total_pos["forward"][PE_size] = 1
							else:
								total_pos["forward"][PE_size] += 1
							
							pos = 0
							for interval in first_Mismatch[:-1]:
								pos += interval + 1
							adj_next_Mismatch = next_Mismatch
							adj_next_Mismatch[0] = next_Mismatch[0] + 80 + gap - pos
							Mismatch = first_Mismatch[:-1] + adj_next_Mismatch
							pos = 0
							for interval in Mismatch[:-1]:
								pos += interval + 1
								mismatch_pos["forward"][PE_size][pos] += 1
				
				elif FLAG == "163":
					first_sequence = Read
					first_MDZ = alignment[12].split(":")[-1]
					first_Mismatch = map(int, re.split(r'[ACGTN]',first_MDZ))
					next_alignment = f.next().strip().split('\t')
					next_sequence = Seq(next_alignment[9])
					next_MDZ = next_alignment[12].split(':')[-1]
					next_Mismatch = map(int, re.split(r'[ACGTN]',next_MDZ))
					next_ID = next_alignment[0]
					next_FLAG = next_alignment[1]
					
					# count reads containing 3' mismatch
					for i in first_Mismatch:
						if i ==0:
							counter += 1
						else:
							break
					reverse_3t_counter[counter] += 1
					sequence = Seq(Read).reverse_complement()
					MisNT = sequence[0-counter:]
					if counter > 0:
						if MisNT not in reverse_3t_NT:
							reverse_3t_NT[MisNT] = 1
						else:
							reverse_3t_NT[MisNT] += 1
					
					# Begin counting mismatches at different positions
					if next_ID == ID and next_FLAG == "83":
						if PE_size <= 160:
							head_length = tail_length = PE_size - 80
							shared_length = 160 - PE_size
							if first_sequence[80 - shared_length:80] == next_sequence[0:shared_length]:
								if PE_size not in mismatch_pos["reverse"]:
									mismatch_pos["reverse"][PE_size] = {}
									for i in range(1, PE_size + 1):
										mismatch_pos["reverse"][PE_size][i] = 0
								if PE_size not in total_pos["reverse"]:
									total_pos["reverse"][PE_size] = 1
								else:
									total_pos["reverse"][PE_size] += 1

								next_Mismatch[0] += head_length
								pos = 0
								adj_next_Mismatch = []
								for interval in next_Mismatch:
									pos += interval + 1
									if pos > 80:
										adj_next_Mismatch.append(interval)
								pos = 0
								for interval in first_Mismatch:
									pre_pos = pos
									pos += interval + 1
									if pos > head_length:
										adj_next_Mismatch[0] -= pre_pos
										break
								Mismatch = first_Mismatch[:-1] + adj_next_Mismatch
								pos = 0
								for interval in Mismatch[::-1][:-1]:
									pos += interval + 1
									mismatch_pos["reverse"][PE_size][pos] += 1
							else:
								continue
						else:
							gap = PE_size - 160
							if PE_size not in mismatch_pos["reverse"]:
								mismatch_pos["reverse"][PE_size] = {}
								for i in range(1, PE_size + 1):
									mismatch_pos["reverse"][PE_size][i] = 0
							if PE_size not in total_pos["reverse"]:
								total_pos["reverse"][PE_size] = 1
							else:
								total_pos["reverse"][PE_size] += 1
							
							pos = 0
							for interval in first_Mismatch[:-1]:
								pos += interval + 1
							adj_next_Mismatch = next_Mismatch
							adj_next_Mismatch[0] = next_Mismatch[0] + 80 + gap - pos
							
							Mismatch = first_Mismatch[:-1] + adj_next_Mismatch
							pos = 0
							for interval in Mismatch[::-1][:-1]:
								pos += interval + 1
								mismatch_pos["reverse"][PE_size][pos] += 1
				else:
					continue
		return (forward_3t_counter, reverse_3t_counter, forward_3t_NT, reverse_3t_NT, mismatch_pos, total_pos)

def update_dict(old_dict, new_dict):
	for i in old_dict:
		if i not in new_dict:
			new_dict[i] = old_dict[i]
		else:
			new_dict[i] += old_dict[i]

def main(arg_input):
	sum_forward_3t_counter = {}
	sum_reverse_3t_counter = {}
	sum_forward_3t_NT = {}
	sum_reverse_3t_NT = {}
	sum_mismatch_pos = {"forward":{}, "reverse":{}}
	sum_total_pos = {"forward":{}, "reverse":{}}
	
	for sam in arg_input:
		mismatch_output = terminal_misincorporation(sam)
		update_dict(mismatch_output[0],sum_forward_3t_counter)
		update_dict(mismatch_output[1],sum_reverse_3t_counter)
		update_dict(mismatch_output[2],sum_forward_3t_NT)
		update_dict(mismatch_output[3],sum_reverse_3t_NT)

		for strand in mismatch_output[4]:
			for size in mismatch_output[4][strand]:
				if size not in sum_mismatch_pos[strand]:
					sum_mismatch_pos[strand][size] = mismatch_output[4][strand][size]
				else:
					for pos in mismatch_output[4][strand][size]:
						if pos not in sum_mismatch_pos[strand][size]:
							sum_mismatch_pos[strand][size][pos] = mismatch_output[4][strand][size][pos]
						else:
							sum_mismatch_pos[strand][size][pos] += mismatch_output[4][strand][size][pos]

		for strand in mismatch_output[5]:
			for size in mismatch_output[5][strand]:
				if size not in sum_total_pos:
					sum_total_pos[strand][size] = mismatch_output[5][strand][size]
				else:
					sum_total_pos[strand][size] += mismatch_output[5][strand][size]

	max_size = max(max(sum_total_pos["forward"], sum_total_pos["reverse"]))

	with open("terminal_misincorporation_count_stat.txt",'w') as f:
		f.write("%s\t%s\t%s\n" % ("Strand","Number of misincorporation","Count"))
		for i in sum_forward_3t_counter:
			f.write("%s\t%d\t%d\n" % ("RDR2 Transcripts", i, sum_forward_3t_counter[i]))
		for i in sum_reverse_3t_counter:
			f.write("%s\t%d\t%d\n" % ("PolIV Transcripts", i, sum_reverse_3t_counter[i]))
	with open("terminal_misincorporation_NT_stat.txt",'w') as f:
		f.write("%s\t%s\t%s\n" % ("Strand","Misincorporated nucleotide","Count"))
		for i in sum_forward_3t_NT:
			f.write("%s\t%s\t%d\n" % ("RDR2 Transcripts", i, sum_forward_3t_NT[i]))
		for i in sum_reverse_3t_NT:
			f.write("%s\t%s\t%d\n" % ("PolIV Transcripts", i, sum_reverse_3t_NT[i]))

	with open("mismatch_position.txt",'w') as f:
		for strand in sum_mismatch_pos:
			for size in sum_mismatch_pos[strand]:
				for pos in sum_mismatch_pos[strand][size]:
					f.write("%s\t%s\t%s\t%d\t%d\t%.4f\n" %
						(strand, size, pos, sum_mismatch_pos[strand][size][pos], sum_total_pos[strand][size], 
							sum_mismatch_pos[strand][size][pos]/sum_total_pos[strand][size]))

	with open("mismatch_position_PolIV",'w') as f:
		f.write("%s\t%s\t%s\n" % ("Strand","Length",'\t'.join(map(str,range(1,max_size)))))				
		for size in sum_mismatch_pos["reverse"]:
			mismatch_count_list = []
			for pos in range(1, size + 1):
				mismatch_count_list.append(sum_mismatch_pos["reverse"][size][pos])
			mismatch_percent_list = map(lambda x: ("%.4f") % (x/sum_total_pos["reverse"][size]), mismatch_count_list)
			f.write("%s\t%s\t%s\n" % ("PolIV Transcripts",size,'\t'.join(mismatch_percent_list)))

	with open("mismatch_position_RDR2",'w') as f:
		f.write("%s\t%s\t%s\n" % ("Strand","Length",'\t'.join(map(str,range(1,max_size)))))
		for size in sum_mismatch_pos["forward"]:
			mismatch_count_list = []
			for pos in range(1, size + 1):
				mismatch_count_list.append(sum_mismatch_pos["forward"][size][pos])
			mismatch_percent_list = map(lambda x: ("%.4f") % (x/sum_total_pos["forward"][size]), mismatch_count_list)
			f.write("%s\t%s\t%s\n" % ("RDR2 Transcripts",size,'\t'.join(mismatch_percent_list)))


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-q','--query',
		help="Query file, SAM format",
		nargs='+',
		dest='input_sam',
		required=True)
	args = parser.parse_args()

	main(args.input_sam)

	
















