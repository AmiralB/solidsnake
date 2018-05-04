### Alignment
configfile: "config.yml"


rule master:
  input:
    expand("{sample}.sam", sample = config["SAMPLES"])


rule alignment:
    input:
        R1 = config["FASTQ_DIR"] + "/" +"{sample}_R1.fastq",
        R2 = config["FASTQ_DIR"] + "/" +"{sample}_R2.fastq"
    output:
        "{sample}.sam"
    shell:
        "bwa mem -R '@RG\\tID:{wildcards.sample}\\tSM:${wildcards.sample}\\tLB:GENOMIQUE\\tPL:ILLUMINA' -t {threads} {config[HG19_PATH]} {input.R1} {input.R2} > {output}"

 