#!/usr/bin/env python
# usage: python PE_trimadapter.py -r1 GSF1696-I-14_S14_R1_001.fastq -r2 GSF1696-I-14_S14_R2_001.fastq -a1 TGGAATTC -a2 GATCGTCG -output_merged 14_trimmed.fq -untrimmed_r1 14_untrimmed_R1.fq -untrimmed_r2 14_untrimmed_R2.fq
# This is a special adapter trimming script for paired end sequencing experiment.
# It is used for the M13 sequencing experiments in the Jasleen et al. paper.
# Adapter sequence for R1 reads: TGGAATTC; Adapter sequence for R2 reads: GATCGTCG

from itertools import izip
from Bio.Seq import Seq
import re
import regex
import argparse


def PE_trimadapter(r1,r2,adapter1,adapter2,trimmed_file,untrimmed_r1,untrimmed_r2):
	trimmed_write = 0
	trimmed_unwrite = 0
	single_trimmed = 0
	untrimmed = 0
	with open(trimmed_file,'w') as f_trimmed, open(untrimmed_r1,'w') as f_untrimmed_r1, open(untrimmed_r2,'w') as f_untrimmed_r2:
		with open(r1,'rb') as f1, open(r2,'rb') as f2:
			for (line1,line2) in izip(f1,f2):
				header1 = line1.strip().split(' ')[0]
				header2 = line2.strip().split(' ')[0]
				read1 = f1.next().strip()
				read2 = f2.next().strip()
				linker1 = f1.next().strip()
				linker2 = f2.next().strip()
				qual1 = f1.next().strip()
				qual2 = f2.next().strip()
				adp1 = regex.findall(r'(.*)(%s){s<=1}' % adapter1, read1)
				adp2 = regex.findall(r'(.*)(%s){s<=1}' % adapter2, read2)
				
				# Short transcripts, in which adapters can be found in both pair-ended reads.
				# So that the paired reads can be trimmed and regard as a single-ended read.
				if adp1 != [] and adp2 != []:
					trimmed_read1 = adp1[0][0]
					trimmed_read2 = adp2[0][0]
					trimmed_qual1 = qual1[0:len(adp1[0][0])]
					trimmed_qual2 = qual2[0:len(adp2[0][0])]
					if len(trimmed_read1) >= 15 and len(trimmed_read2) >= 15 and Seq(trimmed_read1) == Seq(trimmed_read2).reverse_complement():
						f_trimmed.write("%s\n%s\n%s\n%s\n" % (header1, trimmed_read1, linker1, trimmed_qual1))
						trimmed_write += 1
					else:
						trimmed_unwrite += 1
				
				# If adapter sequences (at least 8 bases) cannot be found:
				# Check if the 5' ends of the paired reads perfectly overlap with each other, and the 3' ends contains adapter sequences (0-7 bases).
				# This suggests that the transcript is 73 nt to 80 nt in length. Adapters are trimmed. Paired reads are condensed into single-ended reads.
				# Otherwise, keep the paired reads untrimmed.
				else:
					if adp1 == [] and adp2 == []:
						read1_head = read1[0:72]
						read2_head = read2[0:72]
						overlap1 = re.findall(r'(.*)(%s)(.*)' % Seq(read2_head).reverse_complement(), read1)
						overlap2 = re.findall(r'(.*)(%s)(.*)' % Seq(read1_head).reverse_complement(), read2)
						if overlap1 != [] and overlap2 != []:
							l = len(overlap1[0][-1])
							if overlap1[0][-1] == adapter1[0:l] and overlap2[0][-1] == adapter2[0:l]:
								trimmed_read1 = overlap1[0][0]+overlap1[0][1]
								trimmed_qual1 = qual1[0:len(overlap1[0][0]+overlap1[0][1])]
								f_trimmed.write("%s\n%s\n%s\n%s\n" % (header1, trimmed_read1, linker1, trimmed_qual1))
							else:
								f_untrimmed_r1.write("%s\n%s\n%s\n%s\n" %(header1, read1, linker1, qual1))
								f_untrimmed_r2.write("%s\n%s\n%s\n%s\n" %(header2, read2, linker2, qual2))
						else:
							f_untrimmed_r1.write("%s\n%s\n%s\n%s\n" %(header1, read1, linker1, qual1))
							f_untrimmed_r2.write("%s\n%s\n%s\n%s\n" %(header2, read2, linker2, qual2))
					else:
						f_untrimmed_r1.write("%s\n%s\n%s\n%s\n" %(header1, read1, linker1, qual1))
						f_untrimmed_r2.write("%s\n%s\n%s\n%s\n" %(header2, read2, linker2, qual2))
						single_trimmed += 1

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-r1',
		help="R1 file of paired-end sequencing",
		dest='r1',
		required=True)
	parser.add_argument('-r2',
		help="R2 file of paired-end sequencing",
		dest='r2',
		required=True)
	parser.add_argument('-a1',
		help="sequence of 3' end adapter to be trimmed off. Must be at least 8 nt.",
		dest='a1',
		required=True)
	parser.add_argument('-a2',
		help="reverse complement sequence of 5' end adapter to be trimmed off. Must be at least 8 nt.",
		dest='a2',
		required=True)
	parser.add_argument('-output_merged',
		help="merged reads. These reads are treated as single-ended reads",
		dest='merged',
		required=True)
	parser.add_argument('-untrimmed_r1',
		help="R1 reads that cannot be trimmed",
		dest='untrim_r1',
		required=True)
	parser.add_argument('-untrimmed_r2',
		help="R2 reads that cannot be trimmed",
		dest='untrim_r2',
		required=True)
	args=parser.parse_args()

	PE_trimadapter(args.r1, args.r2, args.a1, args.a2, args.merged, args.untrim_r1, args.untrim_r2)
