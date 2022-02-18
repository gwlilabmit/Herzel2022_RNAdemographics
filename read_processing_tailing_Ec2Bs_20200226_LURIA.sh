#!/bin/sh

GENOME1="/home/herzel/data/annotations/ecoli/E_coli_K12_MG1655_NC_000913.2"
GENOME2="/home/herzel/data/annotations/bsubtilis/bs-chr"

for FASTQ in $(cat $1)
	do
		
		OUT_NAME=$(echo $FASTQ | awk 'BEGIN { FS = "/" } ; { print $NF } ' | sed 's/.fastq//')
		echo $OUT_NAME
		
		# set min length in cutadapt to 1 (--overlap=1)
		# trim A-tail, C-tail, G-tail, T-tail, CCA$
		# CCA (often at tRNA 3' ends) trimming necessary to avoid confusion of single nucleotide A-tails with CCA tails
		cutadapt --cores=0 -a "CCA$" --minimum-length 15 -O 3 --trimmed-only --info-file="$OUT_NAME"_CCA.info -o "$OUT_NAME"_CCA.fastq $FASTQ
		cutadapt --cores=0 -a "CCA$" --minimum-length 15 -O 3 --discard-trimmed -o "$OUT_NAME"_noCCA.fastq $FASTQ
		cutadapt --cores=0 -a "A{38}" --minimum-length 15 -O 1 --trimmed-only --info-file="$OUT_NAME"_A.info -o "$OUT_NAME"_A.fastq "$OUT_NAME"_noCCA.fastq
		cutadapt --cores=0 -a "C{38}" --minimum-length 15 -O 1 --trimmed-only --info-file="$OUT_NAME"_C.info -o "$OUT_NAME"_C.fastq "$OUT_NAME"_noCCA.fastq
		cutadapt --cores=0 -a "G{38}" --minimum-length 15 -O 1 --trimmed-only --info-file="$OUT_NAME"_G.info -o "$OUT_NAME"_G.fastq "$OUT_NAME"_noCCA.fastq
		cutadapt --cores=0 -a "T{38}" --minimum-length 15 -O 1 --trimmed-only --info-file="$OUT_NAME"_T.info -o "$OUT_NAME"_T.fastq "$OUT_NAME"_noCCA.fastq

		cat "$OUT_NAME"_CCA.info | awk 'BEGIN { FS = "\t" } ; $2==0 { print $0 } ' > tmp.info
		mv tmp.info "$OUT_NAME"_CCA.info

		cat "$OUT_NAME"_A.info | awk ' BEGIN { FS = "\t" } ; $2==0 { print $0 } ' > tmp.info
		mv tmp.info "$OUT_NAME"_A.info
		cat "$OUT_NAME"_C.info | awk ' BEGIN { FS = "\t" } ; $2==0 { print $0 } ' > tmp.info
		mv tmp.info "$OUT_NAME"_C.info
		cat "$OUT_NAME"_G.info | awk ' BEGIN { FS = "\t" } ; $2==0 { print $0 } ' > tmp.info
		mv tmp.info "$OUT_NAME"_G.info
		cat "$OUT_NAME"_T.info | awk ' BEGIN { FS = "\t" } ; $2==0 { print $0 } ' > tmp.info
		mv tmp.info "$OUT_NAME"_T.info

		ls "$OUT_NAME"_[ACGT].fastq "$OUT_NAME"_CCA.fastq > inBowtieUn
		
		if [[ "$OUT_NAME" == "Bs"* ]]
		then
			# mapping to B. subtilis
			
			for FILE in $(cat inBowtieUn)
       		do
               		NAME=$(echo $FILE | sed 's/.fastq//')
               		echo $NAME
                	bowtie -v 0 -k 1 --quiet --sam $GENOME2 "$NAME".fastq > "$NAME"_Bs.sam

                	samtools view -S -b "$NAME"_Bs.sam > tmp.bam
                	samtools view -b -F 4 tmp.bam > "$NAME"_Bs.bam

                	rm "$NAME"_Bs.sam

                	bamToBed -i "$NAME"_Bs.bam > "$NAME"_Bs.bed
                	cat "$NAME"_Bs.bed | awk ' { if ($6=="+") print $1"\t"$3-1"\t"$3"\t+" ; if ($6=="-") print $1"\t"$2"\t"$2+1"\t-"} ' | sort | uniq -c | awk ' { print $2"\t"$3"\t"$4"\t"$1"\t"$5 } ' > "$NAME"_Bs.bedgraph
        	       	grep "-" "$NAME"_Bs.bedgraph | cut -f 1-4 > "$NAME"_Bs_r.bedgraph
                	grep "+" "$NAME"_Bs.bedgraph | cut -f 1-4 > "$NAME"_Bs_f.bedgraph

                	rm "$NAME"_Bs.bed "$NAME"_Bs.bedgraph

        	done
		
		elif [[ "$OUT_NAME" == "Ec"* ]]
		then
			# mapping to E. coli
			
			for FILE in $(cat inBowtieUn)
       		do
               		NAME=$(echo $FILE | sed 's/.fastq//')
               		echo $NAME
                	bowtie -v 0 -k 1 --quiet --sam $GENOME1 "$NAME".fastq > "$NAME"_Ec2.sam

                	samtools view -S -b "$NAME"_Ec2.sam > tmp.bam
                	samtools view -b -F 4 tmp.bam > "$NAME"_Ec2.bam

                	rm "$NAME"_Ec2.sam

                	bamToBed -i "$NAME"_Ec2.bam > "$NAME"_Ec2.bed
                	cat "$NAME"_Ec2.bed | awk ' { if ($6=="+") print $1"\t"$3-1"\t"$3"\t+" ; if ($6=="-") print $1"\t"$2"\t"$2+1"\t-"} ' | sort | uniq -c | awk ' { print $2"\t"$3"\t"$4"\t"$1"\t"$5 } ' > "$NAME"_Ec2.bedgraph
        	       	grep "-" "$NAME"_Ec2.bedgraph | cut -f 1-4 > "$NAME"_Ec2_r.bedgraph
                	grep "+" "$NAME"_Ec2.bedgraph | cut -f 1-4 > "$NAME"_Ec2_f.bedgraph

                	rm "$NAME"_Ec2.bed "$NAME"_Ec2.bedgraph

        	done
			
		fi
		
		
		
	
	done

