"""
A snakemake rule file which includes rules for bgzipping, tabix indexing, and merging bcf/vcfs.
Kept in a separate rule file due to maintain visibility of the main, analysis rules in analysis.smk
"""

rule bgzip:
    input:
        "results/bcfs/{sample}.calls.vcf",    
    output:
        "results/bcfs/{sample}.calls.vcf.gz",
    log:
        "logs/bgzip/{sample}.log",
    wrapper:
        "v1.17.4-17-g62b55d45/bio/bgzip"

rule tabix:
    input:
        calls = "results/bcfs/{sample}.calls.vcf.gz",    
    output:
        calls_tbi = "results/bcfs/{sample}.calls.vcf.gz.tbi",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    log:
        "logs/tabix/{sample}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """

rule bcftools_merge:
    input:
        bcfs = expand("results/bcfs/{sample}.calls.vcf.gz", sample=samples),
        idx = expand("results/bcfs/{sample}.calls.vcf.gz.tbi", sample=samples)
    output:
        vcf = "results/vcfs/{dataset}.merged.vcf",
    log:
        "logs/bcftools/merge_{dataset}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    threads: 12
    shell:
        """
        bcftools merge --threads {threads} -o {output.vcf} -O v {input.bcfs} --force-samples 2> {log}
        """

rule bcftools_merge1:
    input:
        bcfs = expand("results/bcfs/{sample}.calls.vcf.gz", sample=samples1),
        idx = expand("results/bcfs/{sample}.calls.vcf.gz.tbi", sample=samples1)
    output:
        vcf = "results/vcfs/{dataset}.1.vcf",
    log:
        "logs/bcftools/merge1_{dataset}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    threads: 12
    shell:
        """
        bcftools merge --threads {threads} -o {output.vcf} -O v {input.bcfs} 2> {log}
        """

rule bcftools_merge2:
    input:
        vcfs = expand("results/bcfs/{sample}.calls.vcf.gz", sample=samples2),
        idx = expand("results/bcfs/{sample}.calls.vcf.gz.tbi", sample=samples2)
    output:
        vcf = "results/vcfs/{dataset}.2.vcf",
    log:
        "logs/bcftools/merge2_{dataset}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    threads: 12
    shell:
        """
        bcftools merge --threads {threads} -o {output.vcf} -O v {input.vcfs} 2> {log}
        """

rule bgzip2:
    input:
        "results/{ref}/vcfs/{dataset}.{n}.vcf",
    output:
        "results/{ref}/vcfs/{dataset}.{n}.vcf.gz",
    log:
        "logs/bgzip/{dataset}.{n}_{ref}.log",
    wrapper:
        "v1.17.4-17-g62b55d45/bio/bgzip"

rule tabix2:
    input:
        vcfgz = "results/vcfs/{dataset}.{n}.vcf.gz",
    output:
        tbi = "results/vcfs/{dataset}.{n}.vcf.gz.tbi"
    log:
        "logs/tabix/{dataset}.{n}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    shell:
        """
        tabix {input.vcfgz} 2> {log}
        """

rule bcftools_merge3:
    input:
        vcf = expand("results/vcfs/{{dataset}}.{n}.vcf.gz", n=[1,2]),
        tbi = expand("results/vcfs/{{dataset}}.{n}.vcf.gz.tbi", n=[1,2]),
    output:
        vcf = "results/vcfs/{dataset}_merged.vcf",
        holder = touch("results/vcfs/.complete.{dataset}.merge_vcfs")
    log:
        "logs/bcftools/merge3_{dataset}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    shell:
        """
        bcftools merge -o {output.vcf} -Ov {input.vcf} 2> {log}
        """
