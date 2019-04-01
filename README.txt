# This is the general analysis pipeline for Singh et al. (2019) paper
# Dependencies: samtools (Version: 1.5), bowtie (version 1.2.2)
# In-house script: PE_trimadapter.py, Transcript_sequences.py, mismatch.py
# GEO accession number: GSE126086
# Please contact Feng Wang (wangfeng@iu.edu) if you have any questions


# Precursor library 1

# Trim adapters, there will be two output files: trimmed reads and untrimmed reads
python PE_trimadapter.py -r1 GSF1696-I-13_S13_R1_001.fastq -r2 GSF1696-I-13_S13_R2_001.fastq -a1 TGGAATTC -a2 GATCGTCG -output_merged 13_trimmed.fq -untrimmed_r1 13_untrimmed_R1.fq -untrimmed_r2 13_untrimmed_R2.fq

# Align trimmed reads to A. thaliana TAIR10 assembly
bowtie -a -v 0 --best --strata Athaliana_167_TAIR10 -q 13_trimmed.fq -S > L13_trimmed_Ath_unmapped_bt1.sam

# Filter out any reads that can be mapped to Arabidopsis genome, convert remaining reads to fq format
samtools view -f 4 L13_trimmed_Ath_unmapped_bt1.sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > L13_trimmed_Ath_unmapped_bt1.fq

# Align remaining reads to M13mp18 reference 'genome'
bowtie -a -v 3 --best --strata M13mp18_Bayou.fa -q L13_trimmed_Ath_unmapped_bt1.fq -S > L13_trimmed_Ath_M13_bt1.sam

# Align untrimmed reads to A. thaliana TAIR10 assembly, convert remaining reads to fq format
bowtie -a -v 0 --best --strata Athaliana_167_TAIR10 -q -1 13_untrimmed_R1.fq -2 13_untrimmed_R2.fq --allow-contain -S > L13_untrimmed_Ath_bt1.sam

# Filter out any reads that can be mapped to Arabidopsis genome
samtools view -f 1 -f 4 -f 8 L13_untrimmed_Ath_bt1.sam |awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > L13_untrimmed_Ath_unmapped_bt1_r1.fq
samtools view -f 1 -f 4 -f 8 L13_untrimmed_Ath_bt1.sam |awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > L13_untrimmed_Ath_unmapped_bt1_r2.fq

# Align remaining untrimmed reads to M13mp18 reference 'genome'
bowtie -a -v 3 --best --strata M13mp18_Bayou.fa -q -1 L13_untrimmed_Ath_unmapped_bt1_r1.fq -2 L13_untrimmed_Ath_unmapped_bt1_r2.fq --allow-contain -S > L13_untrimmed_Ath_M13_bt1.sam

# Assemble RNA transcripts by using the alignment files of both trimmed and untrimmed reads
python Transcript_sequences.py -s L13_trimmed_Ath_M13_bt1.sam -p L13_untrimmed_Ath_M13_bt1.sam -ref M13mp18_Bayou.fa -out_s L13_short.txt -out_p L13_long.txt

# Calculate/count the mismatches in the RNA transcripts
python mismatch_v2_PE.py -q L13_trimmed_Ath_M13_bt1.sam L13_untrimmed_Ath_M13_bt1.sam

# DCL3 treated library 1

# Trim adapters, there will be two output files: trimmed reads and untrimmed reads
python PE_trimadapter.py -r1 GSF1696-I-14_S14_R1_001.fastq -r2 GSF1696-I-14_S14_R2_001.fastq -a1 TGGAATTC -a2 GATCGTCG -output_merged 14_trimmed.fq -untrimmed_r1 14_untrimmed_R1.fq -untrimmed_r2 14_untrimmed_R2.fq

# Align trimmed reads to A. thaliana TAIR10 assembly
bowtie -a -v 0 --best --strata Athaliana_167_TAIR10 -q 14_trimmed.fq -S > L14_trimmed_Ath_unmapped_bt1.sam

# Filter out any reads that can be mapped to Arabidopsis genome, convert remaining reads to fq format
samtools view -f 4 L14_trimmed_Ath_unmapped_bt1.sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > L14_trimmed_Ath_unmapped_bt1.fq

# Align remaining reads to M13mp18 reference 'genome'
bowtie -a -v 3 --best --strata M13mp18_Bayou.fa -q L14_trimmed_Ath_unmapped_bt1.fq -S > L14_trimmed_Ath_M13_bt1.sam

# Align untrimmed reads to A. thaliana TAIR10 assembly, convert remaining reads to fq format
bowtie -a -v 0 --best --strata Athaliana_167_TAIR10 -q -1 14_untrimmed_R1.fq -2 14_untrimmed_R2.fq --allow-contain -S > L14_untrimmed_Ath_bt1.sam

# Filter out any reads that can be mapped to Arabidopsis genome
samtools view -f 1 -f 4 -f 8 L14_untrimmed_Ath_bt1.sam |awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > L14_untrimmed_Ath_unmapped_bt1_r1.fq
samtools view -f 1 -f 4 -f 8 L14_untrimmed_Ath_bt1.sam |awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > L14_untrimmed_Ath_unmapped_bt1_r2.fq

# Align remaining untrimmed reads to M13mp18 reference 'genome'
bowtie -a -v 3 --best --strata M13mp18_Bayou.fa -q -1 L14_untrimmed_Ath_unmapped_bt1_r1.fq -2 L14_untrimmed_Ath_unmapped_bt1_r2.fq --allow-contain -S > L14_untrimmed_Ath_M13_bt1.sam

# Assemble RNA transcripts by using the alignment files of both trimmed and untrimmed reads
python Transcript_sequences.py -s L14_trimmed_Ath_M13_bt1.sam -p L14_untrimmed_Ath_M13_bt1.sam -ref M13mp18_Bayou.fa -out_s L14_short.txt -out_p L14_long.txt

