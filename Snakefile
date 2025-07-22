import pandas as pd
# Import data_table
SAMPLES_INFO=pd.read_csv('data_table.tsv',sep="\t")
# Create new column
SAMPLES_INFO['Sample'] = SAMPLES_INFO["GSM"] + '_' + SAMPLES_INFO["Cell"] + '_' + SAMPLES_INFO["Target"]
SAMPLES_INFO['File'] = SAMPLES_INFO['File'].str.strip()
configfile: "config.yaml"

# Target rule for all rules
rule all:
    input:
        expand("qc/fastqc/{sample}.html", sample=SAMPLES_INFO['Sample']),
        expand("qc/fastqc/{sample}_fastqc.zip", sample=SAMPLES_INFO['Sample']),
        "qc/multiqc/reads.html",
        expand("indexes/{genome}/{genome}.fa.gz", genome = config['genome']),
        expand("indexes/{genome}.{index}.bt2", genome=config["genome"], index=range(1,5)),
        expand("indexes/{genome}.rev.{index}.bt2", genome=config["genome"], index=range(1,3)),
        expand("bams/{sample}_{genome}.bam", sample=SAMPLES_INFO['Sample'], genome=config["genome"]),
        "qc/multiqc/bams.html",
        expand("bams_sorted/{sample}_{genome}.sorted.bam", sample=SAMPLES_INFO['Sample'], genome=config["genome"]),
        expand("bams_sorted/{sample}_{genome}.sorted.bai", sample=SAMPLES_INFO['Sample'], genome=config["genome"]),
        expand("bams_sorted/{sample}_{genome}.sorted.coverage.bw", sample=SAMPLES_INFO['Sample'], genome=config["genome"]),
        multiqc = 'results/multiqc_results.tar.gz',
        bams = 'results/bam_multiqc_results.tar.gz',
        bc ='results/bigwig_results.tar.gz'

#Use wrapper
rule run_fastqc:
    input:
        lambda wildcards: expand("reads/{file_name}",
                                 file_name=SAMPLES_INFO.loc[SAMPLES_INFO['Sample'] == wildcards.sample, 'File'])
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip"
    params:
        ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
           "0.57.0/bio/fastqc"

rule multiqc:
    input:
        expand("qc/fastqc/{sample}.html", sample=SAMPLES_INFO['Sample'])
    output:
        "qc/multiqc/reads.html"
    params:
        ""
    log:
        "logs/multiqc.log"
    wrapper:
        "0.57.0/bio/multiqc"

rule download:
    output:
        expand("indexes/{genome}/{genome}.fa.gz", genome = config['genome'])
    shell:
        "wget http://hgdownload.soe.ucsc.edu/goldenPath/{config[genome]}/bigZips/{config[genome]}.fa.gz -P indexes/{config[genome]}/"

#Was ok
rule bowtie2:
    input:
        expand("indexes/{genome}/{genome}.fa.gz", genome = config["genome"])
    output:
        expand("indexes/{genome}.{index}.bt2", genome=config["genome"], index=range(1,5)),
        expand("indexes/{genome}.rev.{index}.bt2", genome=config["genome"], index=range(1,3))
    params:
        output_name = expand("indexes/{genome}", genome=config["genome"])
    conda:
        "envs/bowtie.yaml"
    shell:
        "bowtie2-build {input} {params.output_name}"

#Wrapper
rule alignment:
    input:
        sample = lambda wildcards: expand("reads/{file_name}",
                                 file_name=SAMPLES_INFO.loc[SAMPLES_INFO['Sample'] == wildcards.sample, 'File'])
    output:
        "bams/{sample}_{genome}.bam"
    log:
        "logs/bowtie2/{sample}_{genome}.log"
    params:
        index = expand("indexes/{genome}", genome=config["genome"]),
        extra=""
    threads: 4
    wrapper:
        "0.58.0/bio/bowtie2/align"

rule multiqc_2:
    input:
        expand("logs/bowtie2/{sample}_{genome}.log", sample=SAMPLES_INFO['Sample'], genome=config["genome"])
    output:
        "qc/multiqc/bams.html"
    shell:
         "multiqc -n {output} -n {input}"

#Wrapper
rule samtools_sort:
    input:
        expand("bams/{sample}_{genome}.bam", sample=SAMPLES_INFO['Sample'], genome=config["genome"])
    output:
        "bams_sorted/{sample}_{genome}.sorted.bam"
    params:
        "-m 4G"
    threads:
        4
    wrapper:
        "0.59.2/bio/samtools/sort"

rule samtools_index:
    input:
        rules.samtools_sort.output
    output:
        "bams_sorted/{sample}_{genome}.sorted.bai"
    params:
        "" # optional params string
    wrapper:
        "0.59.2/bio/samtools/index"

rule coverage:
    input:
         bai = rules.samtools_index.output,
         bam = rules.samtools_sort.output
    output:
          "bams_sorted/{sample}_{genome}.sorted.coverage.bw"
    conda:
          "envs/bamCoverage.yaml"
    shell:
         """
         bamCoverage \
             --bam {input.bam} \
             --outFileName {output} \
             --numberOfProcessors 2
                   """

rule archieve:
    input:
         multiqc = "qc/multiqc/reads.html",
         bams = "qc/multiqc/bams.html",
         bc = expand("bams_sorted/{sample}_{genome}.sorted.coverage.bw", sample=SAMPLES_INFO['Sample'], genome=config["genome"])
    output:
         multiqc = 'results/multiqc_results.tar.gz',
         bams = 'results/bam_multiqc_results.tar.gz',
         bc ='results/bigwig_results.tar.gz'
    shell:
         """
                  tar cvzf {output.multiqc} {input.multiqc} &&
                   tar cvzf {output.bams} {input.bams}
                   tar cvzf {output.bc} {input.bc}
                   """
