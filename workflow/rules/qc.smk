rule index_read_fastqc:
    input:
        sample=rules.bcl_convert.output,
    output:
        index1_qc="results/qc/index-read-qc/I1.html",
        index2_qc="results/qc/index-read-qc/I2.html",
    conda:
        "../envs/AmpSeeker-qc.yaml"
    log:
        "logs/index-read-quality.log"
    shell:
        """
        zcat resources/reads/*I1*.fastq.gz | fastqc stdin --outdir results/qc/index-read-qc/ 2>> {log}
        mv results/qc/index-read-qc/stdin_fastqc.html results/qc/index-read-qc/I1.html 2>> {log}
        mv results/qc/index-read-qc/stdin_fastqc.zip results/qc/index-read-qc/I1.zip 2>> {log}
        
        zcat resources/reads/*I2*.fastq.gz | fastqc stdin --outdir results/index-read-qc/ 2>> {log}
        mv results/qc/index-read-qc/stdin_fastqc.html results/qc/index-read-qc/I2.html 2>> {log}
        mv results/qc/index-read-qc/stdin_fastqc.zip results/qc/index-read-qc/I2.zip 2>> {log}
        """

rule fastp:
    input:
        sample=["resources/reads/{sample}_1.fastq.gz", "resources/reads/{sample}_2.fastq.gz"]
    output:
        trimmed=["results/trimmed-reads/{sample}_1.fastq.gz", "results/trimmed-reads/{sample}_2.fastq.gz"],
        html="results/qc/fastp_reports/{sample}.html",
        json="results/qc/fastp_reports/{sample}.json",
    log:
        "logs/fastp/{sample}.log"
    threads: 4
    wrapper:
        "v1.25.0/bio/fastp"


rule mosdepthCoverage:
  """
  Target per-base coverage with mosdepth
  """
    input:
        bam="results/alignments/{sample}.bam",
        idx="results/alignments/{sample}.bam.bai"
    output:
        "results/coverage/{sample}.per-base.bed.gz",
        "results/coverage/{sample}.mosdepth.summary.txt",
        "results/coverage/{sample}.mosdepth.global.dist.txt",
    log:
        "logs/coverage/{sample}.log"
    threads:4
    conda:
        "../envs/AmpSeeker-cli.yaml"
    params:
        prefix="results/coverage/{sample}",
    shell:
        """
        mosdepth {params.prefix} {input.bam} --fast-mode --threads {threads} 2> {log}
        """

rule bam_stats:
  """
  Calculate mapping statistics with samtools flagstat
  """
    input:
        bam = "results/alignments/{sample}.bam",
        idx = "results/alignments/{sample}.bam.bai"
    output:
        stats = "results/alignments/bamStats/{sample}.flagstat"
    conda:
        "../envs/AmpSeeker-cli.yaml"
    log:
        "logs/BamStats/{sample}.log"
    shell:
        """
        samtools flagstat {input.bam} > {output} 2> {log}
        """

# qualimap analysis for alignment QC 
rule qualimap:
    input:
        bam="results/alignments/{sample}.bam",
    output:
        folder = directory("results/qc/qualimap/{sample}"),
        txt = "results/qc/qualimap/{sample}/genome_results.txt",
    log:
        "logs/qualimap/bamqc/{sample}.log",
    conda:
        "../envs/AmpSeeker-qualimap.yaml"
    shell:
        """
        qualimap bamqc -bam {input.bam} -outdir {output.folder} 2> {log}
        """

rule vcf_stats:
    input:
        vcf = "results/vcfs/targets/{dataset}.merged.vcf"
    output:
        stats = "results/qc/{dataset}.merged.vcf.txt"
    log:
        "logs/vcfStats/{dataset}.log"
    conda:
        "../envs/AmpSeeker-cli.yaml"
    shell:
        """
        bcftools stats {input.vcf} > {output.stats} 2> {log}
        """

rule multiQC:
    input:
        expand("results/qc/fastp_reports/{sample}.json", sample=samples) if config['quality-control']['fastp'] else [],
        expand("results/alignments/bamStats/{sample}.flagstat", sample=samples) if config['quality-control']['stats'] else [],
        expand("results/qc/{dataset}.merged.vcf.txt", dataset=dataset) if config['quality-control']['stats'] else [],
        expand("results/coverage/{sample}.per-base.bed.gz", sample=samples) if config['quality-control']['coverage'] else [],
        expand("results/coverage/{sample}.mosdepth.summary.txt", sample=samples) if config['quality-control']['coverage'] else [],
        expand("results/coverage/{sample}.mosdepth.global.dist.txt", sample=samples) if config['quality-control']['coverage'] else [],
        expand("results/qc/qualimap/{sample}/genome_results.txt", sample=samples) if config['quality-control']['qualimap'] else [],
    output:
        "results/qc/multiqc/multiqc_report.html"
    params:
        extra="--config resources/multiqc.yaml"
    log:
        "logs/multiqc/multiqc.log"
    conda:
        "../envs/AmpSeeker-qc.yaml"
    wrapper: 
        "v2.2.1/bio/multiqc"