configfile: "config.yml"


rule master:
	input:
		expand("{sample}.sam", sample = config["SAMPLES"])


### Alignment

rule alignment:
	input:
		R1 = config["FASTQ_DIR"] + "/" +"{sample}_R1.fastq",
		R2 = config["FASTQ_DIR"] + "/" +"{sample}_R2.fastq"
	output:
		"{sample}.sam"
	log:
		"{sample}.bwa"
	shell:
		'bwa mem -a -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:GENOMIQUE\\tPL:illumina" '
		'-t {threads} {config[HG19_PATH]} {input.R1} {input.R2} > {output} > {log}' 

rule samtobam:
	input:
		"{sample}.sam"
	output:
		"{sample}.bam"
	shell:
		'samtools view -b -@{threads} {input} > {output}'

rule sortbam:
	input:
		"{sample}.bam"
	output:
		"{sample}.sort.bam"
	shell:
		'samtools sort {input} -o {output} -@{threads} -O BAM' #-m 3G améliore le temps de process?


### Prepare Bam
## Mark Duplicates

rule markduplicates:
	input:
		"{sample}.sort.bam"
	output:
		"{sample}.sort.md.bam"
	log:
		metrics = "{sample}.sort.md.metrics",
		stdout  = "{sample}.picard.stdout"

	shell:
		'picard -Xmx4g -XX:ParallelGCThreads={threads} '
		'MarkDuplicates I={input} O={output} METRICS_FILE={log.metrics} '
	 	'VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true 2> {log.stdout} ' #MAX_RECORDS_IN_RAM=5000000 CREATE_MD5_FILE=true READ_NAME_REGEX=null 


##Base Recalibration

rule baserecalibration:
	input:
		"{sample}.sort.md.bam"
	output:
		"{sample}.sort.md.recal.table"
	shell:
		'gatk --java-options "-Xmx4g -XX:ParallelGCThreads={threads}" '
		'BaseRecalibrator -R {config[HG19_PATH]} -I {input} -O {output} '
		'--known-sites {config[DBSNP_PATH]} --known-sites {config[GOLDINDEL_PATH]} '
		'--use-original-qualities --read-validation-stringency LENIENT '

rule applyBQSR:
	input:
		bam = "{sample}.sort.md.bam",
		table = "{sample}.sort.md.recal.table"
	output:
		"{sample}.sort.md.recal.bam"
	shell:
		'gatk --java-options "-Xmx4g -XX:ParallelGCThreads={threads}" '
		'ApplyBQSR -R {config[HG19_PATH]} -I {input.bam} -O {output} -bqsr {input.table} '
		'--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 '
		'--use-original-qualities --read-validation-stringency LENIENT ' #--add-output-sam-program-record --create-output-bam-md5


### Variant Calling
## Haplotypecaller

rule haplotypecaller:
	input:
		"{sample}.sort.md.recal.bam"
	output:
		"{sample}.haplotypecaller.snp_indel.g.vcf"
	shell:
		'gatk --java-options "-Xmx4g -XX:ParallelGCThreads={threads}" '
		'HaplotypeCaller -R {config[HG19_PATH]} -I {input} -O {output} -L {config[REFSEQ_PATH]} -D {config[DBSNP_PATH]} '
		'-ip 100 -ERC GVCF ' #-stand-call-conf 10.0 -nct/-nt --min-pruning 3 --dont-use-soft-clipped-bases

rule combineGVCFs:
	input:
		expand("{sample}.haplotypecaller.snp_indel.g.vcf", sample = config["SAMPLES"])
	output:
		"all.haplotypecaller.snp_indel.g.vcf"
	shell:
		'gatk --java-options "-Xmx4g -XX:ParallelGCThreads={threads}" '
		'CombineGVCFs -R {config[HG19_PATH]} -V {input} -O {output} -L {config[REFSEQ_PATH]} -D {config[DBSNP_PATH]} '

rule genotypeGVCFs:
	input:
		"all.haplotypecaller.snp_indel.g.vcf"
	output:
		"all.haplotypecaller.snp_indel.vcf"
	shell:
		'gatk --java-options "-Xmx4g -XX:ParallelGCThreads={threads}" '
		'GenotypeGVCFs -R {config[HG19_PATH]} -V {input} -O {output} -L {config[REFSEQ_PATH]} -D {config[DBSNP_PATH]} '
		'-ip 100 -stand-call-conf 10.0 ' #-nct/-nt


## Freebayes

# rule fixbai:
# 	input:
# 		"{sample}.sort.md.bam"
# 	output:
# 		"{sample}.sort.md.fx.bam"
# 	shell:
# 		'{config[FIXBAI_PATH]} -i {input} > {output} '

rule freebayes:
	input:
		expand("{sample}.sort.md.bam", sample = config["SAMPLES"]) #expand("{sample}.sort.md.fx.bam", sample = config["SAMPLES"])
	output:
		"all.freebayes.snp_indel.vcf"
	shell:
		'freebayes --fasta-reference {config[HG19_PATH]} -v {output} -b {input} -t {config[REFSEQ_PATH]} '
		'--min-alternate-count 2 --min-alternate-qsum 40 --pvar 0.0001 --use-mapping-quality --posterior-integration-limits 1,3 '
		'--genotype-variant-threshold 4 --site-selection-max-iterations 3 --genotyping-max-iterations 25 --max-complex-gap 3 '

# rule sortvcf:
# 	input:
# 		"all.freebayes.snp_indel.vcf"
# 	output:
# 		"all.freebayes.snp_indel.sort.vcf"
# 	shell:
# 		'{config[VCFSORT_PATH]} {config[HG19_DICT_PATH]} {input} > {output} '

# rule sortvcf:
# 	input:
# 		"all.freebayes.snp_indel.vcf"
# 	output:
# 		"all.freebayes.snp_indel.sort.vcf"
# 	shell:
# 		'picard -Xmx4g -XX:ParallelGCThreads={threads} '
#		'SortVcf R={config[HG19_PATH]} I={input} O={output} SD={confi[HG19_DICT_PATH]} '


##Hard Filter

rule variantfiltration:
	input:
		"all.{caller}.snp_indel.vcf"
	output:
		"all.{caller}.snp_indel.filtered.vcf"
	shell:
		'gatk --java-options "-Xmx4g -XX:ParallelGCThreads={threads}" ' 
		'VariantFiltration -R {config[HG19_PATH]} -V {input} -O {output} --clusterWindowSize 10 '
		'--filter-expression "MQ0 >= 4 && ((MQ0 / (1.0*DP)) > 0.1)" --filter-name "HARD_TO_VALIDATE" '
		'--filter-expression "DP < 5" --filter-name "LOW_COVERAGE" '
		'--filter-expression "QUAL < 30.0" --filter-name "VERY_LOW_QUAL" '
		'--filter-expression "QUAL > 30.0 && QUAL < 50.0" --filter-name "LOW_QUAL" '
		'--filter-expression "QD < 1.5" --filter-name "LOW_QD" '
		'--filterExpression "FS > 60.0 " --filterName "FisherStrandBias" ' #voir si erreur (tp)

rule selectSNP:
	input:
		"all.{caller}.snp_indel.filtered.vcf"
	output:
		"all.{caller}.snp.vcf"
	shell:
		'gatk --java-options "-Xmx4g -XX:ParallelGCThreads={threads}" '
		'SelectVariants -R {config[HG19_PATH]} -V {input} -O {output} -select-type SNP '

rule selectINDEL:
	input:
		"all.{caller}.snp_indel.filtered.vcf"
	output:
		"all.{caller}.indel.vcf"
	shell:
		'gatk --java-options "-Xmx4g -XX:ParallelGCThreads={threads}" '
		'SelectVariants -R {config[HG19_PATH]} -V {input} -O {output} -select-type INDEL '
		
rule SNPfiltration:
	input:
		"all.{caller}.snp.vcf"
	output:
		"all.{caller}.snp.filtered.vcf"
	shell:
		'gatk --java-options "-Xmx4g -XX:ParallelGCThreads={threads}" '
		'VariantFiltration -R {config[HG19_PATH]} -V {input} -O {output} --clusterWindowSize 10 '
		'--filter-expression "DP < 10 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "SNP_FILTER" '

rule INDELfiltration:
	input:
		"all.{caller}.indel.vcf"
	output:
		"all.{caller}.indel.filtered.vcf"
	shell:
		'gatk --java-options "-Xmx4g -XX:ParallelGCThreads={threads}" '
		'VariantFiltration -R {config[HG19_PATH]} -V {input} -O {output} --clusterWindowSize 10 '
		'--filter-expression "DP < 10 || QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "INDEL_FILTER" '

rule mergefiltration:
	input:
		snp = "all.{caller}.snp.filtered.vcf",
		indel = "all.{caller}.indel.filtered.vcf"
	output:
		"all.{caller}.snp_indel.filtered2.vcf"
	shell:
		'picard -Xmx4g -XX:ParallelGCThreads={threads} '
		'MergeVcfs R={config[HG19_PATH]} I={input.snp} I={input.indel} O={output} '


### Prepare Annotation
## Combine

rule combineVCFs:
	input:
		haplo = "all.haplotypecaller.snp_indel.filtered2.vcf",
		fb = "all.freebayes.snp_indel.filtered2.vcf"
	output:
		"all.combined.snp_indel.vcf"
	shell:
		'picard -Xmx4g -XX:ParallelGCThreads={threads} '
		'MergeVcfs R={config[HG19_PATH]} I={input.haplo} I={input.fb} O={output} '


## Clean

rule VCFnormalization: #decompose?
	input:
		"all.combined.snp_indel.vcf"
	output:
		"all.combined.snp_indel.clean.vcf"
	shell:
		'vt normalize {input} -r {config[HG19_PATH]} | vt uniq -o {output} '


## Relatedness

rule relatedness:
	input:
		"all.combined.snp_indel.clean.vcf"
	output:
		"vcftools_relatedness2"
	shell:
		"vcftools --vcf {input} --relatedness2 --out {output}"


### Annotation

rule annotation:
	input:
		"all.combined.snp_indel.clean.vcf"
	output:
		"all.ann.vcf"
	shell:
		'snpEff -Xmx4g -XX:ParallelGCThreads={threads} '
		'{config[DBSNPEFF]} {input} > {output} '