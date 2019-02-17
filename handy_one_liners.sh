#!/bin/bash -l

#This command takes sequences in fasta files and removes line breaks
awk '!/^>/ { printf "%s", $0; n = "\n"} /^>/ { print n $0; n="" } END { printf "%s", n } ' barley_transcripts.fa > barley_transcripts_final.fa
#Then this command is used with the fasta file formatted with the above command to pull out transcript sequences by ID
grep -A 1 -f 7DO_SAMvTcAXM_DE_Genes.txt barley_transcripts_final.fa > 7DO_SAMvTcAXM_DE_Genes.fa

#Use to get position for a specific gene ID
grep 'HORVU3HG0397100' Hv_IBSC_PGSB_v2_Dez2015_final.gtfgre

#Use to get positions for a list of gene IDs
grep -f DE_Genes_Seedling_SAMvAXM.txt Hv_IBSC_PGSB_v2_Dez2015_final.gtf > Pos_DE_Genes_Seedling_SAMvAXM.txt

#Use to get BAC IDs if it falls in appropriate range - this will get BACs for all chromosomes
awk '$3 <= 108142516 && $4 >= 108142516' 150831_barley_pseudomolecules_AGP_file.tsv

#Use to get fasta sequence for a specific BAC
grep -A 1 'HVVMRXALLHA0555P07_C5' 150831_pseudomolecule_sequence_contigs.fasta > HVVMRXALLHA0555P07_C5.fa

#Use to get fasta sequence for a specific transcript 
grep -A 1 'HORVU3HG0115000' barley_transcripts_final.fa > HORVU3HG0115000.fa

#The positions 181309065 and 337279253 are the start and stop of the cul2 introgression (area that contains the markers that are in all recombinants at least). Introgressed region may be longer than this.
#Used this command to get a list of genes in the cul2 introgression. To get start and stop positions of the each gene, I first made a copy of the output. For the original output, I used Excel to sort alphabetically by #gene ID and remove duplicates. Then I sorted the copied output first alphabetically by gene ID and then end position largest to smallest. Then I got rid of duplicates in this list. I intersected the two resulting lists #and kept the start position from the original output and the end position from the copied output.
awk '$4 >= 181309065 && $5 <= 337279253' Hv_IBSC_PGSB_v2_Dez2015_final.gtf > genes_in_cul2_region.txt

#Used to get annotations for all genes in the cul2 introgression.
grep -f gene_IDs_in_cul2_region.txt Hv_v1_AHRD_nonTEgenes__15-12-2015-15-54_finalids.csv > Annotations_genes_in_cul2_region.txt

#Used to get positions of genes in the cul2 introgression
grep -f gene_IDs_in_cul2_region.txt Hv_IBSC_PGSB_v2_Dez2015_final.gtf > Pos_Genes_in_cul2_region.txt

#Used to pull BACs that contain genes in cul2 introgression
awk 'FNR==NR{start[NR]=$1+0; stop[NR]=$2+0;next} {for (i in start)if ($3<=start[i] && $4>=stop[i]){print;next}}' cul2_gene_pos.txt 150831_barley_pseudomolecules_AGP_file.tsv > info.txt

#Pull BAC IDs for all SNPs in a file
awk 'FNR==NR{start[NR]=$4+0; stop[NR]=$4+0;next} {for (i in start)if ($3<=start[i] && $4>=stop[i]){print;next}}' 50K_GBS_Merged_Filtered_and_Imputed_TAGS_Chr1.hmp.txt Chr1_AGP_file.tsv > Chr1_BACS.txt

#Unix commands used to pull info out of alignment summaries from Tophat
echo ${x}
grep -m2 'Input' ${x} | tail -n1
grep 'Aligned pairs:' ${x}
grep -m3 'have multiple alignments' ${x} | tail -n1
grep 'discordant alignments' ${x}

#ran in a .sh script called grep.sh (run by sh grep.sh), and just copy pasted values into spreadsheet
for x in *summary.txt;
do
echo ${x}
done

#convert VCF and GTF files to sorted BED files using BEDOPS
~/software/bin/vcf2bed < HapMap_GBS_50K_Merged_TAGS.vcf > sorted_GBS_50K_Merged_TAGS.vcf.bed
~/software/bin/gtf2bed < Hv_IBSC_PGSB_v2_Dez2015_final.gtf > Hv_IBSC_PGSB_v2_Dez2015_final.gtf.bed

#change formatting of chromosome numbers to numeric
sed -i 's/chr1H/1/g' Hv_IBSC_PGSB_v2_Dez2015_final.gtf.bed

#change order of columns and get rid of column #5
awk '{print $8,$6,$7,$9,$1,$2,$3,$4}' file1.txt > file2.txt

