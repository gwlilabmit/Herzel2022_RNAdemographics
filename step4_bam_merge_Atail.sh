#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 4                    
#SBATCH --mail-type=END
#SBATCH --mail-user=herzel@mit.edu
#SBATCH -o /home/herzel/log_LURIA/merge_Ec2A_bam."%j".out
#SBATCH -e /home/herzel/log_LURIA/merge_Ec2A_bam."%j".err
#SBATCH --export=ALL
hostname
#############################################

module add samtools/1.5
module add bedtools/2.25.0

echo Today date is:
date
echo Your current directory is:

cd /home/herzel/data/data/3end_seq/210903Li_poly_cond/tailing/
ls  /home/herzel/data/data/3end_seq/210903Li_poly_cond/bam/*_trim_noScBs_Ec2.bam > inBAM_SI

# merge bam files after A-tail read mapping and without

for BAM in $(cat inBAM_SI )
	do
		OUT_NAME=$(echo $BAM | awk 'BEGIN { FS = "/" } ; { print $NF } ' | sed 's/_trim_noScBs_Ec2.bam//')		
		echo $OUT_NAME
		
		BAM_NT=$(echo $BAM | sed 's/bam/tailing\/bam/' | sed 's/_trim_noScBs_Ec2.bam/_trim_unScEc2BsOSI_A_Ec2.bam/')
		echo $BAM_NT

		samtools merge "$OUT_NAME"_trim_Ec2_allA.bam $BAM $BAM_NT
		# convert from bam to bedgraph
		bamToBed -i "$OUT_NAME"_trim_Ec2_allA.bam > "$OUT_NAME"_trim_Ec2_allA.bed
		cat "$OUT_NAME"_trim_Ec2_allA.bed | awk ' { if ($6=="+") print $1"\t"$3-1"\t"$3"\t+" ; if ($6=="-") print $1"\t"$2"\t"$2+1"\t-"} ' | sort | uniq -c | awk ' { print $2"\t"$3"\t"$4"\t"$1"\t"$5 } ' > "$OUT_NAME"_trim_Ec2_allA.bedgraph
		cat "$OUT_NAME"_trim_Ec2_allA.bedgraph | awk ' $5=="-" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_trim_Ec2_allA_r.bedgraph
		cat "$OUT_NAME"_trim_Ec2_allA.bedgraph | awk ' $5=="+" { print $0 } ' | cut -f 1-4 > "$OUT_NAME"_trim_Ec2_allA_f.bedgraph
		
		rm "$OUT_NAME"_trim_Ec2_allA.bed "$OUT_NAME"_trim_Ec2_allA.bedgraph

	done

