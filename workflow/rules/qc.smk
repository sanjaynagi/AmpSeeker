rule index_read_fastqc:
    input:
        sample=rules.bcl_convert.output,
    output:
        index1_qc="results/index-read-qc/I1.html",
        index2_qc="results/index-read-qc/I2.html",
    conda:
        "../envs/AmpSeeker-qc.yaml"
    log:
        "logs/index-read-quality.log"
    shell:
        """
        zcat resources/reads/*I1*.fastq.gz | fastqc stdin --outdir results/index-read-qc/ 2>> {log}
        mv results/index-read-qc/stdin_fastqc.html results/index-read-qc/I1.html 2>> {log}
        mv results/index-read-qc/stdin_fastqc.zip results/index-read-qc/I1.zip 2>> {log}
        
        zcat resources/reads/*I2*.fastq.gz | fastqc stdin --outdir results/index-read-qc/ 2>> {log}
        mv results/index-read-qc/stdin_fastqc.html results/index-read-qc/I2.html 2>> {log}
        mv results/index-read-qc/stdin_fastqc.zip results/index-read-qc/I2.zip 2>> {log}
        """

rule fastp:
    input:
        sample=["resources/reads/{sample}_1.fastq.gz", "resources/reads/{sample}_2.fastq.gz"]
    output:
        trimmed=["results/trimmed-reads/{sample}_1.fastq.gz", "results/trimmed-reads/{sample}_2.fastq.gz"],
        html="results/fastp_reports/{sample}.html",
        json="results/fastp_reports/{sample}.json",
        logs="logs/fastp/{sample}.log"
    log:
        "logs/fastp/{sample}.log"
    threads: 4
    wrapper:
        "v1.25.0/bio/fastp"

rule multiQC:
    input:
        expand("results/fastp_reports/{sample}.json", sample=samples),
        expand("results/alignments/bamStats/{sample}.flagstat", sample=samples),
        expand("results/coverage/{sample}.per-base.bed.gz", sample=samples),
        expand("results/coverage/{sample}.mosdepth.summary.txt", sample=samples),
        expand("results/coverage/{sample}.mosdepth.global.dist.txt", sample=samples),
        expand("results/qualimap/{sample}/genome_results.txt", sample=samples),
        expand("results/vcfs/stats/{dataset}.merged.vcf.txt", dataset=dataset),
    output:
        "results/multiqc/multiqc_report.html"
    params:
        extra="--config config/multiqc.yaml"
    log:
        "logs/multiqc/multiqc.log"
    conda:
        "../envs/AmpSeeker-qc.yaml"
    wrapper: 
        "v2.2.1/bio/multiqc"


rule mosdepthCoverage:
  """
  Target per-base coverage with mosdepth
  """
    input:
        bam="results/alignments/{sample}.bam",
        idx="results/alignments/{sample}.bam.bai"
    output:
        "results/coverage/{sample}.per-base.bed.gz"
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
        folder = directory("results/qualimap/{sample}"),
        txt = "results/qualimap/{sample}/genome_results.txt",
    log:
        "logs/qualimap/bamqc/{sample}.log",
    conda:
        "../envs/AmpSeeker-qualimap.yaml"
    shell:
        """
        qualimap bamqc -bam {input.bam} -outdir {output.folder} -outformat PDF:HTML 2> {log}
        """

rule vcf_stats:
    input:
        vcf = "results/vcfs/{dataset}.merged.vcf"
    output:
        stats = "results/vcfs/stats/{dataset}.merged.vcf.txt"
    log:
        "logs/vcfStats/{dataset}.log"
    conda:
        "../envs/AmpSeeker-cli.yaml"
    shell:
        """
        bcftools stats {input.vcf} > {output.stats} 2> {log}
        """
