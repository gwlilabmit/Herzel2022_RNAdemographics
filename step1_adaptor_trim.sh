#!/bin/bash
#SBATCH -N 1                     
#SBATCH -n 4                    
#SBATCH --mail-type=END
#SBATCH --mail-user=herzel@mit.edu
#SBATCH -o /home/herzel/log_LURIA/adaptor_trim."%j".out
#SBATCH -e /home/herzel/log_LURIA/adaptor_trim."%j".err
#SBATCH --export=ALL
hostname
#############################################

module add bowtie/1.2
module add python/2.7.13

echo Today date is:
date

# this is a representative script for processing individual 3' end sequencing datasets
# here two datasets are processed: pre- and post- Kasugamycin treatment (Figure 5D)
cd /home/herzel/data/data/3end_seq/210903Li_poly_cond/

cp ~/rawdata/210903Li/210903Li_D21-8276*.fastq Ec_WT_preK_3.fastq
cp ~/rawdata/210903Li/210903Li_D21-8277*.fastq Ec_WT_postK_3.fastq

python /home/herzel/data/scripts/2_removeLinker_09052018_linkerRemovalOnly.py Ec_WT_preK_3.fastq
python /home/herzel/data/scripts/2_removeLinker_09052018_linkerRemovalOnly.py Ec_WT_postK_3.fastq

gzip ~/rawdata/210903Li/210903Li_D21-8276*.fastq
gzip ~/rawdata/210903Li/210903Li_D21-8277*.fastq

mv Ec_WT_preK_3.f-stripped.fq Ec_WT_preK_3_trim.fastq
mv Ec_WT_postK_3.f-stripped.fq Ec_WT_postK_3_trim.fastq
