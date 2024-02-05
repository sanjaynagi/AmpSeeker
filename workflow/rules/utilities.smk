"""
A snakemake rule file which includes rules for bgzipping, tabix indexing, and merging bcf/vcfs.
Kept in a separate rule file due to maintain visibility of the main, analysis rules in analysis.smk
"""


rule bgzip:
    input:
        vcf="results/vcfs/{call_type}/{sample}.calls.vcf",
    output:
        vcfgz="results/vcfs/{call_type}/{sample}.calls.vcf.gz",
    log:
        log="logs/bgzip/{call_type}/{sample}.log",
    wrapper:
        "v1.17.4-17-g62b55d45/bio/bgzip"


rule tabix:
    input:
        calls="results/vcfs/{call_type}/{sample}.calls.vcf.gz",
    output:
        calls_tbi="results/vcfs/{call_type}/{sample}.calls.vcf.gz.tbi",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    log:
        "logs/tabix/{call_type}/{sample}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """


rule bcftools_merge:
    input:
        vcfs=expand("results/vcfs/{{call_type}}/{sample}.calls.vcf.gz", sample=samples),
        idx=expand(
            "results/vcfs/{{call_type}}/{sample}.calls.vcf.gz.tbi", sample=samples
        ),
    output:
        vcf="results/vcfs/{call_type}/{dataset}.merged.vcf",
    log:
        "logs/bcftools/{call_type}/merge_{dataset}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    threads: 12
    shell:
        """
        bcftools merge --threads {threads} -o {output.vcf} -O v {input.vcfs} --force-samples 2> {log}
        """


rule bcftools_merge1:
    input:
        vcfs=expand("results/vcfs/{call_type}/{sample}.calls.vcf.gz", sample=samples1),
        idx=expand(
            "results/vcfs/{call_type}/{sample}.calls.vcf.gz.tbi", sample=samples1
        ),
    output:
        vcf="results/vcfs/{call_type}/{dataset}.1.vcf",
    log:
        "logs/bcftools/{call_type}/merge1_{dataset}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    threads: 12
    shell:
        """
        bcftools merge --threads {threads} -o {output.vcf} -O v {input.vcfs} 2> {log}
        """


rule bcftools_merge2:
    input:
        vcfs=expand("results/vcfs/{call_type}/{sample}.calls.vcf.gz", sample=samples2),
        idx=expand(
            "results/vcfs/{call_type}/{sample}.calls.vcf.gz.tbi", sample=samples2
        ),
    output:
        vcf="results/vcfs/{call_type}/{dataset}.2.vcf",
    log:
        "logs/bcftools/{call_type}/merge2_{dataset}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    threads: 12
    shell:
        """
        bcftools merge --threads {threads} -o {output.vcf} -O v {input.vcfs} 2> {log}
        """


rule bgzip2:
    input:
        vcf="results/vcfs/{call_type}/{dataset}.{n}.vcf",
    output:
        vcfgz="results/vcfs/{call_type}/{dataset}.{n, [0-9]}.vcf.gz",
    log:
        "logs/bgzip/{call_type}/{dataset}.{n}.log",
    wrapper:
        "v1.17.4-17-g62b55d45/bio/bgzip"


rule tabix2:
    input:
        vcfgz="results/vcfs/{call_type}/{dataset}.{n}.vcf.gz",
    output:
        tbi="results/vcfs/{call_type}/{dataset}.{n, [/d]}.vcf.gz.tbi",
    log:
        "logs/tabix/{call_type}/{dataset}.{n}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    shell:
        """
        tabix {input.vcfgz} 2> {log}
        """


rule bcftools_merge3:
    input:
        vcf=expand("results/vcfs/{{call_type}}/{{dataset}}.{n}.vcf.gz", n=[1, 2]),
        tbi=expand("results/vcfs/{{call_type}}/{{dataset}}.{n}.vcf.gz.tbi", n=[1, 2]),
    output:
        vcf="results/vcfs/{call_type}/{dataset}.merged.vcf",
        holder=touch("results/vcfs/{call_type}/.complete.{dataset}.merge_vcfs"),
    log:
        "logs/bcftools/{call_type}/merge3_{dataset}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    shell:
        """
        bcftools merge -o {output.vcf} -Ov {input.vcf} 2> {log}
        """
