#! /bin/bash
#$ -S /bin/bash
#$ -M lemoinehugo1@gmail.com
#$ -N snakefile
#$ -m bea

DIR="result_test3"




#Lancement environnement virtuel
source /local/env/envconda3.sh
source activate solidsnake
cd project/solidsnake/

#Lancement pipeline
snakemake --cluster "/usr/local/sge/bin/lx-amd64/qsub -pe make {threads} -V" -p --jobs 64 --latency-wait 120 all.combined.ann.vcf all.combined

#Rangement
mkdir $DIR
mkdir $DIR/log && mv *.log *.metrics $DIR/log/
mkdir $DIR/bam && mv *.sam *.bam *.table *.bai $DIR/bam/
mkdir -p $DIR/vcf/haplotypecaller && mv *.haplotypecaller.* $DIR/vcf/haplotypecaller/
mkdir -p $DIR/vcf/freebayes && mv *.freebayes.* $DIR/vcf/freebayes/
mkdir -p $DIR/vcf/combined && mv *.combined.* $DIR/vcf/combined/
mkdir $DIR/annotation && mv snpEff_* $DIR/annotation/

#Multiqc
cd $DIR
multiqc . -o ./multiqc_report.html