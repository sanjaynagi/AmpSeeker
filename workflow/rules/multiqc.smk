rule FastQC:
    input:
        reads = expand(["results/reads/{sample}_1.fastq.gz","results/reads/{sample}_2.fastq.gz"], sample=samples),
    output:
        html = expand(["results/fastqc/{sample}_1_fastqc.html", "results/fastqc/{sample}_2_fastqc.html"], sample=samples),
        zip = expand(["results/fastqc/{sample}_1_fastqc.zip", "results/fastqc/{sample}_2_fastqc.zip"], sample=samples)
    log:
        "logs/fastqc/fastq.log"
    conda:
        "../envs/AmpSeeker-qc.yaml"
    threads: 4
    params:
        outdir="--outdir results/fastqc/",
    shell:
        """
        fastqc {input} {params.outdir} -t {threads} 2> {log}
        """
    # wrapper:
    #     "v1.25.0/bio/fastqc"

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
        expand(["results/fastqc/{sample}_1_fastqc.html","results/fastqc/{sample}_2_fastqc.html",
                "results/fastqc/{sample}_1_fastqc.zip", "results/fastqc/{sample}_2_fastqc.zip",
                "results/alignments/bamStats/{sample}.flagstat",
                "results/qualimap/{sample}",
                #"results/qualimap/{sample}.genome_fraction_coverage.txt",
                #"results/qualimap/{sample}.mapped_reads_gc-content_distribution.txt",
                #"results/qualimap/{sample}.genome_results.txt",
                #"results/qualimap/{sample}.coverage_histogram.txt",
                "results/coverage/{sample}.per-base.bed.gz",
                "results/wholegenome/coverage/windowed/{sample}.regions.bed.gz"], sample=samples),
        expand("results/vcfs/stats/{dataset}.merged.vcf.txt", dataset=dataset)
    output:
        "results/multiqc/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        "v1.25.0/bio/multiqc"
