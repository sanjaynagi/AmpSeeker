"""
A snakemake rule file which includes rules for bgzipping, tabix indexing, and merging bcf/vcfs.
Kept in a separate rule file due to maintain visibility of the main, analysis rules in analysis.smk
"""

rule bgzip:
    input:
        "results/{ref}/bcfs/{sample}.calls.vcf",    
    output:
        "results/{ref}/bcfs/{sample}.calls.vcf.gz",
    log:
        "logs/bgzip/{sample}_{ref}.log",
    wrapper:
        "v1.17.4-17-g62b55d45/bio/bgzip"

rule tabix:
    input:
        calls = "results/{ref}/bcfs/{sample}.calls.vcf.gz",    
    output:
        calls_tbi = "results/{ref}/bcfs/{sample}.calls.vcf.gz.tbi",
    conda:
        "../envs/AmpSeq.yaml"
    log:
        "logs/tabix/{sample}_{ref}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """

# wont work with over 1000 files 
rule bcftools_merge1:
    input:
        bcfs = expand("results/{{ref}}/bcfs/{sample}.calls.vcf.gz", sample=samples1),
        idx = expand("results/{{ref}}/bcfs/{sample}.calls.vcf.gz.tbi", sample=samples1)
    output:
        vcf = "results/{ref}/vcfs/AgamDaoLSTM1.vcf",
    log:
        "logs/bcftools_merge_{ref}.log",
    conda:
        "../envs/AmpSeq.yaml"
    threads: 12
    shell:
        """
        bcftools merge --threads {threads} -o {output.vcf} -O v {input.bcfs} 2> {log}
        """

rule bcftools_merge2:
    input:
        vcfs = expand("results/{{ref}}/bcfs/{sample}.calls.vcf.gz", sample=samples2),
        idx = expand("results/{{ref}}/bcfs/{sample}.calls.vcf.gz.tbi", sample=samples2)
    output:
        vcf = "results/{ref}/vcfs/AgamDaoLSTM2.vcf",
    log:
        "logs/bcftools_merge2_{ref}.log",
    conda:
        "../envs/AmpSeq.yaml"
    threads: 12
    shell:
        """
        bcftools merge --threads {threads} -o {output.vcf} -O v {input.vcfs} 2> {log}
        """


rule bgzip2:
    input:
        "results/{ref}/vcfs/AgamDaoLSTM{n}.vcf",
    output:
        "results/{ref}/vcfs/AgamDaoLSTM{n}.vcf.gz",
    log:
        "logs/bgzip/main_{n}_{ref}.log",
    wrapper:
        "v1.17.4-17-g62b55d45/bio/bgzip"

rule tabix2:
    input:
        vcfgz = "results/{ref}/vcfs/AgamDaoLSTM{n}.vcf.gz",
    output:
        tbi = "results/{ref}/vcfs/AgamDaoLSTM{n}.vcf.gz.tbi"
    log:
        "logs/tabix/main_{n}_{ref}.log",
    conda:
        "../envs/AmpSeq.yaml"
    shell:
        """
        tabix {input.vcfgz} 2> {log}
        """

rule bcftools_merge3:
    input:
        vcf = expand("results/{{ref}}/vcfs/AgamDaoLSTM{n}.vcf.gz", n=[1,2]),
        tbi = expand("results/{{ref}}/vcfs/AgamDaoLSTM{n}.vcf.gz.tbi", n=[1,2])
    output:
        vcf = "results/{ref}/vcfs/AgamDaoLSTM_merged.vcf",
    log:
        "logs/bcftools_merge3_{ref}.log",
    conda:
        "../envs/AmpSeq.yaml"
    shell:
        """
        bcftools merge -o {output.vcf} -Ov {input.vcf} 2> {log}
        """
