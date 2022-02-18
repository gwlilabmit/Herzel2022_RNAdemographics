#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 4                    
#SBATCH --mail-type=END
#SBATCH --mail-user=herzel@mit.edu
#SBATCH -o /home/herzel/log_LURIA/mapping_tail."%j".out
#SBATCH -e /home/herzel/log_LURIA/mapping_tail."%j".err
#SBATCH --export=ALL
hostname
#############################################

module add bowtie/1.2
module add samtools/1.5
module add bedtools/2.25.0
module add python/2.7.13
module add cutadapt/1.16

echo Today date is:
date
echo Your current directory is:

# created tailing folder before running script
# created fastq_unmapped folder and moved all unmapped files there before running script
cd /home/herzel/data/data/3end_seq/210903Li_poly_cond/tailing/

ls  /home/herzel/data/data/3end_seq/210903Li_poly_cond/fastq_unmapped/*_unScEc2BsOSI.fastq > inFASTQ

sh ~/scripts_LURIA/read_processing_tailing_Ec2Bs_20200226_LURIA.sh inFASTQ
