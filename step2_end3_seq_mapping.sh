#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 4                    
#SBATCH --mail-type=END
#SBATCH --mail-user=herzel@mit.edu
#SBATCH -o /home/herzel/log_LURIA/mapping."%j".out
#SBATCH -e /home/herzel/log_LURIA/mapping."%j".err
#SBATCH --export=ALL
hostname
#############################################

module add bowtie/1.2
module add samtools/1.5
module add bedtools/2.25.0
module add fastxtoolkit/0.0.13

echo Today date is:
date
echo Your current directory is:

cd /home/herzel/data/data/3end_seq/210903Li_poly_cond/
pwd

ls fastq_trim/*_trim.fastq > inLS

# mapping with bowtie

GENOME1="/home/herzel/data/annotations/ecoli/E_coli_K12_MG1655_NC_000913.2"
GENOME2="/home/herzel/data/annotations/bsubtilis/bs-chr"
GENOME3="/home/herzel/data/annotations/o199/o199"
GENOME4="/home/herzel/data/annotations/IVT_spike_in/IVT_spike_in"
GENOME5="/home/herzel/data/annotations/Scer3/Scer3"

for FILE in $(cat inLS)
	do
		OUT_NAME=$(echo $FILE | awk 'BEGIN { FS = "/" } ; { print $NF } ' | sed 's/.fastq//')
		fastx_trimmer -Q 33 -f 2 -i $FILE -o "$OUT_NAME"_t1.fastq	
		awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) > 14) {print header, seq, qheader, qseq}}' < "$OUT_NAME"_t1.fastq > tmp.fastq
		mv tmp.fastq "$OUT_NAME"_t1_g14.fastq
		rm "$OUT_NAME"_t1.fastq		

		# map to all possible targets 
		bowtie -v 0 -k 1 --quiet --sam $GENOME5 "$OUT_NAME"_t1_g14.fastq  --un "$OUT_NAME"_unSc.fastq > "$OUT_NAME"_Sc.sam

		bowtie -v 0 -k 1 --quiet --sam $GENOME2 "$OUT_NAME"_unSc.fastq --un "$OUT_NAME"_unScBs.fastq > tmp.sam
		bowtie -v 0 -k 1 --quiet --sam $GENOME1 "$OUT_NAME"_unScBs.fastq --un "$OUT_NAME"_unScBsEc2.fastq > "$OUT_NAME"_noScBs_Ec2.sam
		bowtie -v 0 -k 1 --quiet --sam $GENOME3 "$OUT_NAME"_unScBsEc2.fastq --un "$OUT_NAME"_unScEc2BsO.fastq > "$OUT_NAME"_o199.sam
		bowtie -v 0 -k 1 --quiet --sam $GENOME4 "$OUT_NAME"_unScEc2BsO.fastq --un "$OUT_NAME"_unScEc2BsOSI.fastq > "$OUT_NAME"_SI.sam 
		
		bowtie -v 0 -k 1 --quiet --sam $GENOME1 "$OUT_NAME"_unSc.fastq --un "$OUT_NAME"_unScEc2.fastq > tmp.sam
		bowtie -v 0 -k 1 --quiet --sam $GENOME2 "$OUT_NAME"_unScEc2.fastq --un "$OUT_NAME"_unScEc2Bs.fastq > "$OUT_NAME"_noScEc2_Bs.sam

		rm tmp.sam

		# convert sam to bam and remove any unmapped reads to reduce file size	
		samtools view -S -b "$OUT_NAME"_Sc.sam > tmp.bam
		samtools view -b -F 4 tmp.bam > "$OUT_NAME"_Sc.bam
		samtools view -S -b "$OUT_NAME"_noScBs_Ec2.sam > tmp.bam
		samtools view -b -F 4 tmp.bam > "$OUT_NAME"_noScBs_Ec2.bam
		samtools view -S -b "$OUT_NAME"_o199.sam > tmp.bam
		samtools view -b -F 4 tmp.bam > "$OUT_NAME"_o199.bam
		samtools view -S -b "$OUT_NAME"_SI.sam > tmp.bam
		samtools view -b -F 4 tmp.bam > "$OUT_NAME"_SI.bam
		samtools view -S -b "$OUT_NAME"_noScEc2_Bs.sam > tmp.bam
		samtools view -b -F 4 tmp.bam > "$OUT_NAME"_noScEc2_Bs.bam
		
		rm *.sam
		
		# convert to bedgraph to get single nt resolution
		bamToBed -i "$OUT_NAME"_noScEc2_Bs.bam > "$OUT_NAME"_noScEc2_Bs.bed
		cat "$OUT_NAME"_noScEc2_Bs.bed | awk ' { if ($6=="+") print $1"\t"$3-1"\t"$3"\t+" ; if ($6=="-") print $1"\t"$2"\t"$2+1"\t-"} ' | sort | uniq -c | awk ' { print $2"\t"$3"\t"$4"\t"$1"\t"$5 } ' > "$OUT_NAME"_noScEc2_Bs.bedgraph
		cat "$OUT_NAME"_noScEc2_Bs.bedgraph | awk ' $5=="-" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_noScEc2_Bs_r.bedgraph
		cat "$OUT_NAME"_noScEc2_Bs.bedgraph | awk ' $5=="+" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_noScEc2_Bs_f.bedgraph
		
		bamToBed -i "$OUT_NAME"_noScBs_Ec2.bam > "$OUT_NAME"_noScBs_Ec2.bed
		cat "$OUT_NAME"_noScBs_Ec2.bed | awk ' { if ($6=="+") print $1"\t"$3-1"\t"$3"\t+" ; if ($6=="-") print $1"\t"$2"\t"$2+1"\t-"} ' | sort | uniq -c | awk ' { print $2"\t"$3"\t"$4"\t"$1"\t"$5 } ' > "$OUT_NAME"_noScBs_Ec2.bedgraph
		cat "$OUT_NAME"_noScBs_Ec2.bedgraph | awk ' $5=="-" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_noScBs_Ec2_r.bedgraph
		cat "$OUT_NAME"_noScBs_Ec2.bedgraph | awk ' $5=="+" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_noScBs_Ec2_f.bedgraph
		
		bamToBed -i "$OUT_NAME"_o199.bam > "$OUT_NAME"_o199.bed
		cat "$OUT_NAME"_o199.bed | awk ' { if ($6=="+") print $1"\t"$3-1"\t"$3"\t+" ; if ($6=="-") print $1"\t"$2"\t"$2+1"\t-"} ' | sort | uniq -c | awk ' { print $2"\t"$3"\t"$4"\t"$1"\t"$5 } ' > "$OUT_NAME"_o199.bedgraph
		cat "$OUT_NAME"_o199.bedgraph | awk ' $5=="-" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_o199_r.bedgraph
		cat "$OUT_NAME"_o199.bedgraph | awk ' $5=="+" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_o199_f.bedgraph
		
		bamToBed -i "$OUT_NAME"_SI.bam > "$OUT_NAME"_SI.bed
		cat "$OUT_NAME"_SI.bed | awk ' { if ($6=="+") print $1"\t"$3-1"\t"$3"\t+" ; if ($6=="-") print $1"\t"$2"\t"$2+1"\t-"} ' | sort | uniq -c | awk ' { print $2"\t"$3"\t"$4"\t"$1"\t"$5 } ' > "$OUT_NAME"_SI.bedgraph
		cat "$OUT_NAME"_SI.bedgraph | awk ' $5=="-" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_SI_r.bedgraph
		cat "$OUT_NAME"_SI.bedgraph | awk ' $5=="+" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_SI_f.bedgraph
		
		bamToBed -i "$OUT_NAME"_Sc.bam > "$OUT_NAME"_Sc.bed
        cat "$OUT_NAME"_Sc.bed | awk ' { if ($6=="+") print $1"\t"$3-1"\t"$3"\t+" ; if ($6=="-") print $1"\t"$2"\t"$2+1"\t-"} ' | sort | uniq -c | awk ' { print $2"\t"$3"\t"$4"\t"$1"\t"$5 } ' > "$OUT_NAME"_Sc.bedgraph
        cat "$OUT_NAME"_Sc.bedgraph | awk ' $5=="-" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_Sc_r.bedgraph
        cat "$OUT_NAME"_Sc.bedgraph | awk ' $5=="+" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_Sc_f.bedgraph

		rm *.bed "$OUT_NAME"_pLac-CUA.bedgraph "$OUT_NAME"_SI.bedgraph   "$OUT_NAME"_o199.bedgraph "$OUT_NAME"_Sc.bedgraph "$OUT_NAME"_noScBs_Ec2.bedgraph "$OUT_NAME"_noScEc2_Bs.bedgraph
		rm "$OUT_NAME"_unSc.fastq "$OUT_NAME"_unScBs.fastq "$OUT_NAME"_unScBsEc2.fastq "$OUT_NAME"_unScEc2BsO.fastq "$OUT_NAME"_unScEc2.fastq "$OUT_NAME"_unScEc2Bs.fastq "$OUT_NAME"_unpLac-CUA.fastq
		rm tmp.bam
		gzip "$OUT_NAME"_t1_g14.fastq

	done
