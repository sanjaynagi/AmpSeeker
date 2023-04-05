rule FastQC:
    input:
        reads = expand(["results/reads/{sample}_1.fastq.gz","results/reads/{sample}_2.fastq.gz"], sample=samples),
    output:
        html = expand(["results/fastqc/{sample}_1_fastqc.html", "results/fastqc/{sample}_2_fastqc.html"], sample=samples),
        zip = expand(["results/fastqc/{sample}_1_fastqc.zip", "results/fastqc/{sample}_2_fastqc.zip"], sample=samples)
    log:
        "logs/fastqc/fastq.log"
    conda:
        "../envs/AmpSeeker-cli-lock.yaml"
    threads: 4
    params:
        outdir="--outdir results/fastqc/",
    shell:
        """
        fastqc {input} {params.outdir} -t {threads} 2> {log}
        """

rule vcfStats:
    input:
        vcf = "results/vcfs/{dataset}.merged.vcf"
    output:
        stats = "results/vcfs/stats/{dataset}.merged.vcf.txt"
    log:
        "logs/vcfStats/{dataset}.log"
    conda:
        "../envs/AmpSeeker-cli-lock.yaml"
    shell:
        """
        bcftools stats {input.vcf} > {output.stats} 2> {log}
        """

rule multiQC:
    input:
        expand(["results/fastqc/{sample}_1_fastqc.html","results/fastqc/{sample}_2_fastqc.html"], sample=samples),
        expand(["results/fastqc/{sample}_1_fastqc.zip", "results/fastqc/{sample}_2_fastqc.zip"], sample=samples),
        expand("results/alignments/bamStats/{sample}.flagstat", sample=samples),
        expand("results/coverage/{sample}.per-base.bed.gz", sample=samples),
        expand("results/wholegenome/coverage/windowed/{sample}.regions.bed.gz", sample=samples),
        expand("results/vcfs/stats/{dataset}.merged.vcf.txt", dataset=dataset)
    output:
        "results/multiqc/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    conda:
        "../envs/AmpSeeker-cli-lock.yaml"
    shell:
        """
        multiqc results -o results/multiqc/ -f {input} 2> {log}
        """
