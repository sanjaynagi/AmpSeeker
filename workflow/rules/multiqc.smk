rule fastp:
    input:
        sample=["results/reads/{sample}_1.fastq.gz", "results/reads/{sample}_2.fastq.gz"]
    output:
        trimmed=["results/reads/trimmed/{sample}_1.fastq.gz", "results/reads/trimmed/{sample}_2.fastq.gz"],
        html="results/fastp_reports/{sample}.html",
        json="results/fastp_reports/{sample}.json",
    log:
        "logs/fastp/{sample}.log"
    threads: 4
    wrapper:
        "v1.25.0/bio/fastp"

rule vcfStats:
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

rule multiQC:
    input:
        expand("results/fastp_reports/{sample}.html", sample=samples),
        expand("results/fastp_reports/{sample}.json", sample=samples),
        expand("results/alignments/bamStats/{sample}.flagstat", sample=samples),
        expand("results/coverage/{sample}.per-base.bed.gz", sample=samples),
        expand("results/wholegenome/coverage/windowed/{sample}.regions.bed.gz", sample=samples),
        expand("results/vcfs/stats/{dataset}.merged.vcf.txt", dataset=dataset)
    output:
        "results/multiqc/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        "v1.25.0/bio/multiqc"
