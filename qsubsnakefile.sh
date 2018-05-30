#! /bin/bash
# Les commentraires qui commencent par '#$' sont
# interpretes par SGE comme des options en ligne

# Shell a utiliser pour l'execution du job
#$ -S /bin/bash

#$ -o /ngs/datagen/diag-genet/MEX/TEST/QXT_MEX_Fam2618HPE/RESULTS_GATK4.log

#$ -cwd

#$ -N snakefile

# Utilisateur a avertir
#$ -M hugo.lemoine@etudiant.univ-rennes1.fr

# Avertir au debut (b)egin, a la fin (e)nd, a l'eliminaton (a)bort et a la suspension (s)uspend d'un job
#$ -m bea

# Export de toutes les variables d'environnement
#$ -V

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