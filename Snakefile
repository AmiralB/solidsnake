configfile: 'config.yml'

rule master:
	input:
		expand('{sample}.sam', sample = config['SAMPLES']),
		'all.combined.ann.vcf',
		'all.relatedness2'


### Alignment

# rule bcltofastq:


# # rule fastqc:
# 	input:
# 		config['FASTQ_DIR'] + '/' + '{sample}_R1.fastq'
#  	output:
#  		'fastQC/{afile}_fastqc.zip'
#  	shell:
#  		'fastqc -o fastQC {input}'

rule unzip:
	input:
		config['FASTQ_DIR'] + '/' + '{filename}.fastq.gz'
	output:
		config['FASTQ_DIR'] + '/' + '{filename}.fastq'
	log:
		'{filename}.unzip.log'
	#benchmark:
	#	'.decompression.benchmark.txt'
	shell:
		'gunzip {input} > {output} '
		#'2> {log} '

'''
Alignement d'une paire de fastq en sam
'''
rule alignment:
	input:
		R1 = config['FASTQ_DIR'] + '/' + '{sample}_R1.fastq',# gz
		R2 = config['FASTQ_DIR'] + '/' + '{sample}_R2.fastq'
	output:
		'{sample}.sam'
	log:
		'{sample}.bwa.log'	
	#benchmark:
	#	'{sample}.bwa.benchmark.txt'
	threads:
		12
	message:
		"Alignment of fastq"
	shell:
		'bwa mem -a -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:GENOMIQUE\\tPL:illumina" '
		'-t {threads} {config[HG19_PATH]} {input.R1} {input.R2} > {output} '
		#'2> {log} '


rule sortbam:
	input:
		'{sample}.sam'
	output:
		'{sample}.sort.bam'
	#benchmark:
	#	'{sample}.sortbam.benchmark.txt'
	threads:
		12
	shell:
		'samtools sort {input} -o {output} -@{threads} -O BAM ' 


### Prepare Bam
## Mark Duplicates

rule markduplicates:
	input:
		'{sample}.sort.bam'
	output:
		'{sample}.sort.md.bam'
	log:
		metrics = '{sample}.sort.md.metrics',
		errout  = '{sample}.md.log'
	#benchmark:
	#	'{sample}.md.benchmark.txt'
	threads:
		12	
	shell:
		'picard {config[JAVA_ARGS]} -XX:ParallelGCThreads={threads} '
		'MarkDuplicates I={input} O={output} METRICS_FILE={log.metrics} '
		'VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true '
		#'2> {log.errout} '  


##Base Recalibration

rule baserecalibration:
	input:
		'{sample}.sort.md.bam'
	output:
		'{sample}.sort.md.recal.table'
	log:
		'{sample}.br.log'
	#benchmark:
	#	'{sample}.br.benchmark.txt'
	threads:
		12
	shell:
		'gatk --java-options "{config[JAVA_ARGS]} -XX:ParallelGCThreads={threads}" '
		'BaseRecalibrator -R {config[HG19_PATH]} -I {input} -O {output} '
		'--known-sites {config[DBSNP_PATH]} --known-sites {config[GOLDINDEL_PATH]} '
		'--use-original-qualities --read-validation-stringency LENIENT '
		#'2> {log} '


rule applyBQSR:
	input:
		bam = '{sample}.sort.md.bam',
		table = '{sample}.sort.md.recal.table'
	output:
		'{sample}.sort.md.recal.bam'
	log:
		'{sample}.bqsr.log'
	#benchmark:
	#	'{sample}.bqsr.benchmark.txt'
	threads:
		12
	shell:
		'gatk --java-options "{config[JAVA_ARGS]} -XX:ParallelGCThreads={threads}" '
		'ApplyBQSR -R {config[HG19_PATH]} -I {input.bam} -O {output} -bqsr {input.table} '
		'--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 '
		'--use-original-qualities --read-validation-stringency LENIENT ' 
		#'2> {log} '


### Variant Calling
## Haplotypecaller

rule haplotypecaller:
	input:
		'{sample}.sort.md.recal.bam'
	output:
		'{sample}.haplotypecaller.snp_indel.g.vcf'
	log:
		'{sample}.haplotypecaller.log'
	#benchmark:
	#	'{sample}.haplotypecaller.benchmark.txt'
	threads:
		12
	shell:
		'gatk --java-options "{config[JAVA_ARGS]} -XX:ParallelGCThreads={threads}" '
		'HaplotypeCaller -R {config[HG19_PATH]} -I {input} -O {output} -L {config[REFSEQ_PATH]} -D {config[DBSNP_PATH]} '
		'-ip 100 -ERC GVCF '
		#'2> {log} '


rule combineGVCFs:
	input:
		expand('{sample}.haplotypecaller.snp_indel.g.vcf', sample = config['SAMPLES'])
	output:
		'all.haplotypecaller.snp_indel.g.vcf'
	log:
		'all.combineGVCFs.log'
	#benchmark:
	#	'all.cGVCF.benchmark.txt'
	threads:
		12
	shell:
		'gatk --java-options "{config[JAVA_ARGS]} -XX:ParallelGCThreads={threads}" '
		'CombineGVCFs -R {config[HG19_PATH]} -V {input[0]} -V {input[1]} -O {output} -L {config[REFSEQ_PATH]} -D {config[DBSNP_PATH]} '
		#'2> {log} '


rule genotypeGVCFs:
	input:
		'all.haplotypecaller.snp_indel.g.vcf'
	output:
		'all.haplotypecaller.snp_indel.vcf'
	log:
		'all.genotypeGVCFs.log'
	#benchmark:
	#	'all.gGVCF.benchmark.txt'
	threads:
		12
	shell:
		'gatk --java-options "{config[JAVA_ARGS]} -XX:ParallelGCThreads={threads}" '
		'GenotypeGVCFs -R {config[HG19_PATH]} -V {input} -O {output} -L {config[REFSEQ_PATH]} -D {config[DBSNP_PATH]} '
		'-ip 100 -stand-call-conf 10.0 '
		#'2> {log} ' 


## Freebayes

rule freebayes:
	input:
		expand('{sample}.sort.md.bam', sample = config['SAMPLES'])
	output:
		'all.freebayes.snp_indel.unsorted.vcf'
	#benchmark:
	#	'all.freebayes.benchmark.txt'
	shell:
		'freebayes --fasta-reference {config[HG19_PATH]} -v {output} -b {input} -t {config[REFSEQ_PATH]} '
		'--min-alternate-count 2 --min-alternate-qsum 40 --pvar 0.0001 --use-mapping-quality --posterior-integration-limits 1,3 '
		'--genotype-variant-threshold 4 --site-selection-max-iterations 3 --genotyping-max-iterations 25 --max-complex-gap 3 '


rule sortvcf:
	input:
		'all.freebayes.snp_indel.unsorted.vcf'
	output:
		'all.freebayes.snp_indel.vcf'
	log:
		'all.fbsortvcf.log'
	#benchmark:
	#	'all.fbsortvcf.benchmark.txt'
	threads:
		12
	shell:
		'picard {config[JAVA_ARGS]} -XX:ParallelGCThreads={threads} '
		'SortVcf R={config[HG19_PATH]} I={input} O={output} SD={config[HG19_DICT_PATH]} '
		#'2> {log} '


##Hard Filter

rule variantfiltration:
	input:
		'all.{caller}.snp_indel.vcf'
	output:
		'all.{caller}.snp_indel.filtered.vcf'
	log:
		'all.{caller}.variantfiltration.log'
	#benchmark:
	#	'all.variantfilt.benchmark.txt'
	threads:
		12	
	params : 
		 filter_args = " ".join(["--filter-expression '{}' --filter-name '{}'".format(v,i) for i,v in config["FILTER"].items()])
	shell:
		'gatk --java-options "{config[JAVA_ARGS]} -XX:ParallelGCThreads={threads}" ' 
		'VariantFiltration -R {config[HG19_PATH]} -V {input} -O {output} --cluster-window-size 10 '
		'{params.filter_args} '
		#'2> {log} '


rule selectSNP:
	input:
		'all.{caller}.snp_indel.filtered.vcf'
	output:
		'all.{caller}.snp.vcf'
	log:
		'all.{caller}.selectSNP.log'
	#benchmark:
	#	'all.selectSNP.benchmark.txt'
	threads:
		12	
	shell:
		'gatk --java-options "{config[JAVA_ARGS]} -XX:ParallelGCThreads={threads}" '
		'SelectVariants -R {config[HG19_PATH]} -V {input} -O {output} -select-type SNP '
		#'2> {log} '


rule selectINDEL:
	input:
		'all.{caller}.snp_indel.filtered.vcf'
	output:
		'all.{caller}.indel.vcf'
	log:
		'all.{caller}.selectINDEL.log'
	#benchmark:
	#	'all.selectINDEL.benchmark.txt'
	threads:
		12	
	shell:
		'gatk --java-options "{config[JAVA_ARGS]} -XX:ParallelGCThreads={threads}" '
		'SelectVariants -R {config[HG19_PATH]} -V {input} -O {output} -select-type INDEL '
		#'2> {log} '


rule SNPfiltration:
	input:
		'all.{caller}.snp.vcf'
	output:
		'all.{caller}.snp.filtered.vcf'
	log:
		'all.{caller}.SNPfiltration.log'
	#benchmark:
	#	'all.SNPfilt.benchmark.txt'
	threads:
		12	
	shell:
		'gatk --java-options "{config[JAVA_ARGS]} -XX:ParallelGCThreads={threads}" '
		'VariantFiltration -R {config[HG19_PATH]} -V {input} -O {output} --cluster-window-size 10 ' #TODO
		'--filter-expression "DP < 10 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "SNP_FILTER" '
		#'2> {log} '


rule INDELfiltration:
	input:
		'all.{caller}.indel.vcf'
	output:
		'all.{caller}.indel.filtered.vcf'
	log:
		'all.{caller}.INDELfiltration.log'
	#benchmark:
	#	'all.INDELfilt.benchmark.txt'
	threads:
		12	
	shell:
		'gatk --java-options "{config[JAVA_ARGS]} -XX:ParallelGCThreads={threads}" '
		'VariantFiltration -R {config[HG19_PATH]} -V {input} -O {output} --cluster-window-size 10 ' #TODO
		'--filter-expression "DP < 10 || QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "INDEL_FILTER" '
		#'2> {log} '


rule mergefiltration:
	input:
		snp = 'all.{caller}.snp.filtered.vcf',
		indel = 'all.{caller}.indel.filtered.vcf'
	output:
		'all.{caller}.snp_indel.final.vcf'
	log:
		'all.{caller}.mergefiltration.log'
	#benchmark:
	#	'all.mergefilt.benchmark.txt'
	threads:
		12	
	shell:
		'picard {config[JAVA_ARGS]} -XX:ParallelGCThreads={threads} '
		'MergeVcfs R={config[HG19_PATH]} I={input.snp} I={input.indel} O={output} '
		#'2> {log} '


### Prepare Annotation
## Combine

rule combineVCFs:
	input:
		haplo = 'all.haplotypecaller.snp_indel.final.vcf',
		fb    = 'all.freebayes.snp_indel.final.vcf'
	output:
		'all.snp_indel.vcf'
	log:
		'all.cVCF.log'
	#benchmark:
	#	'all.cVCF.benchmark.txt'
	threads:
		12
	shell:
		'picard {config[JAVA_ARGS]} -XX:ParallelGCThreads={threads} '
		'MergeVcfs R={config[HG19_PATH]} I={input.haplo} I={input.fb} O={output} '
		#'2> {log} '


## Clean

rule allelicdecomposition:
	input:
		'all.snp_indel.vcf'
	output:
		'all.snp_indel.decomposed.vcf'
	log:
		'all.decomposition.log'
	#benchmark:
	#	'all.decomposition.benchmark.txt'
	shell:
		'vt decompose {input} -o {output} '
		#'2> {log} '




rule VCFnormalization:
	input:
		'all.snp_indel.vcf'
	output:
		'all.snp_indel.norm.vcf'
	log:
		'all.VCFnormalization.log'
	#benchmark:
	#	'all.VCFnorm.benchmark.txt'
	shell:
		'vt normalize {input} -r {config[HG19_PATH]} | vt uniq - -o {output} '
		#'2>> {log} '
# 2> {log}


### Annotation

rule annotation:
	input:
		'all.snp_indel.norm.vcf'
	output:
		'all.snp_indel.norm.ann.vcf'
	#benchmark:
	#	'all.snpeff.benchmark.txt'
	threads:
		12
	shell:
		'snpEff {config[JAVA_ARGS]} -XX:ParallelGCThreads={threads} '
		'hg19 {input} -t > {output} '
# add other annotation databases (vcf) if neccessary


### Post-Process

rule relatedness:
	input:
		'{prefix}.snp_indel.norm.vcf'
	output:
		'{prefix}.relatedness2'
	log:
		'{prefix}.relatedness2.log'
	#benchmark:
	#	'{prefix}.relatedness2.benchmark.txt'
	shell:
		'vcftools --vcf {input} --relatedness2 --out {wildcards.prefix} '
		#'2> {log} '



### Conversion

rule bamtocram:
	input:
		'{filename}.bam'
	output:
		'{filename}.cram'
	log:
		'{filename}.bamtocram.log'
	#benchmark:
	#	'{filename}.bamtocram.benchmark.txt'
	shell:
		'samtools view -C -T {config[HG19_PATH]} {input} -o {output} '
		#'2> {log} '


# rule compression:
# 	input:
# 		'{filename}'
# 	output:
# 		'{filename}.tgt.gz'
# 	log:
# 		'.compression.log'
# 	#benchmark:
# 	#	'.compression.benchmark.txt'
# 	shell:
# 		'gzip {input} > {output} '
# 		#'2> {log} '


# rule decompression:
# 	input:
# 		'{filename}.gz'
# 	output:
# 		'.'
# 	log:
# 		'.decompression.log'
# 	#benchmark:
# 	#	'.decompression.benchmark.txt'
# 	shell:
# 		'gunzip {input} > {output} '
# 		#'2> {log} '