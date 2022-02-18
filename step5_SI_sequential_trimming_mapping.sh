#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 4                    
#SBATCH --mail-type=END
#SBATCH --mail-user=herzel@mit.edu
#SBATCH -o /home/herzel/log_LURIA/SI_tail."%j".out
#SBATCH -e /home/herzel/log_LURIA/SI_tail."%j".err
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

# this script deals with trimming and mapping of reads with non-templated nucleotides at 3' ends of reference RNAs

# created SI_trim folder before running script
# created fastq_unmapped folder and moved all unmapped files there before running script
cd /home/herzel/data/data/3end_seq/210903Li_poly_cond/SI_trim/

ls  /home/herzel/data/data/3end_seq/210903Li_poly_cond/fastq_unmapped/*_unScEc2BsOSI.fastq > inFASTQ

GENOME4="/home/herzel/data/annotations/IVT_spike_in/IVT_spike_in"

for FASTQ in $(cat inFASTQ )
	do
		OUT_NAME=$(echo $FASTQ | awk 'BEGIN { FS = "/" } ; { print $NF } ' | sed 's/.fastq//')		
		echo $OUT_NAME
		
		i="0"
		while [ "$i" -lt "5" ]
			do
				i=$(expr $i + 1)
				echo $i
				# trim with fastx_trimmer by 1-5 nucleotides
				# keep unmapped reads & trim next nucleotide etc
				if [ $i -eq 1 ]
				then
			#		echo "i equals 1"
					fastx_trimmer -Q 33 -t 1 -m 15 -i $FASTQ -o "$OUT_NAME"_t1.fastq
				else
			#		echo "i does not equal 1"
					fastx_trimmer -Q 33 -t 1 -m 15 -i "$OUT_NAME"_t"$(expr $i - 1)"unSI.fastq -o "$OUT_NAME"_t"$i".fastq
				fi
				
				# after each iteration map to spike-in reference
				bowtie -v 0 -k 1 --quiet --sam $GENOME4 "$OUT_NAME"_t"$i".fastq --un "$OUT_NAME"_t"$i"unSI.fastq > "$OUT_NAME"_t"$i"_SI.sam
				samtools view -S -b "$OUT_NAME"_t"$i"_SI.sam > tmp.bam
				samtools view -b -F 4 tmp.bam > "$OUT_NAME"_t"$i"_SI.bam

			done
	
		# at end merge bam files
		samtools merge "$OUT_NAME"_t1-5_SI.bam "$OUT_NAME"_t1_SI.bam "$OUT_NAME"_t2_SI.bam "$OUT_NAME"_t3_SI.bam "$OUT_NAME"_t4_SI.bam "$OUT_NAME"_t5_SI.bam
		
		# convert from bam to bedgraph
		bamToBed -i "$OUT_NAME"_t1-5_SI.bam > "$OUT_NAME"_t1-5_SI.bed
		cat "$OUT_NAME"_t1-5_SI.bed | awk ' { if ($6=="+") print $1"\t"$3-1"\t"$3"\t+" ; if ($6=="-") print $1"\t"$2"\t"$2+1"\t-"} ' | sort | uniq -c | awk ' { print $2"\t"$3"\t"$4"\t"$1"\t"$5 } ' > "$OUT_NAME"_t1-5_SI.bedgraph
		cat "$OUT_NAME"_t1-5_SI.bedgraph | awk ' $5=="-" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_t1-5_SI_r.bedgraph
		cat "$OUT_NAME"_t1-5_SI.bedgraph | awk ' $5=="+" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_t1-5_SI_f.bedgraph
		
		rm "$OUT_NAME"_t1-5_SI.bed "$OUT_NAME"_t1-5_SI.bedgraph tmp.bam *_SI.sam *_t[1-5].fastq
	
	done



