#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 4                    
#SBATCH --mail-type=END
#SBATCH --mail-user=herzel@mit.edu
#SBATCH -o /home/herzel/log_LURIA/merge_SI_bam."%j".out
#SBATCH -e /home/herzel/log_LURIA/merge_SI_bam."%j".err
#SBATCH --export=ALL
hostname
#############################################

module add samtools/1.5
module add bedtools/2.25.0

echo Today date is:
date
echo Your current directory is:

cd /home/herzel/data/data/3end_seq/210903Li_poly_cond/SI_trim/bam/
ls  /home/herzel/data/data/3end_seq/210903Li_poly_cond/bam/*_trim_SI.bam > inBAM_SI

# merge bam files after spike-in/ reference RNA mapping with non-templated nucleotides and without

for BAM in $(cat inBAM_SI )
	do
		OUT_NAME=$(echo $BAM | awk 'BEGIN { FS = "/" } ; { print $NF } ' | sed 's/_trim_SI.bam//')		
		echo $OUT_NAME
		
		BAM_NT=$(echo $BAM | sed 's/bam/SI_trim\/bam/' | sed 's/_trim_SI.bam/_trim_unScEc2BsOSI_t1-5_SI.bam/')
		echo $BAM_NT

		samtools merge "$OUT_NAME"_trim_SI_all.bam $BAM $BAM_NT
		
		# convert from bam to bedgraph
		bamToBed -i "$OUT_NAME"_trim_SI_all.bam > "$OUT_NAME"_trim_SI_all.bed
		cat "$OUT_NAME"_trim_SI_all.bed | awk ' { if ($6=="+") print $1"\t"$3-1"\t"$3"\t+" ; if ($6=="-") print $1"\t"$2"\t"$2+1"\t-"} ' | sort | uniq -c | awk ' { print $2"\t"$3"\t"$4"\t"$1"\t"$5 } ' > "$OUT_NAME"_trim_SI_all.bedgraph
		cat "$OUT_NAME"_trim_SI_all.bedgraph | awk ' $5=="-" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_trim_SI_all_r.bedgraph
		cat "$OUT_NAME"_trim_SI_all.bedgraph | awk ' $5=="+" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_trim_SI_all_f.bedgraph
		
		rm "$OUT_NAME"_trim_SI_all.bed "$OUT_NAME"_trim_SI_all.bedgraph

	done

