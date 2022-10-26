"""
A snakemake rule file which includes rules for bgzipping, tabix indexing, and merging bcf/vcfs.
Kept in a separate rule file due to maintain visibility of the main, analysis rules in analysis.smk
"""

rule bgzip:
    input:
        calls = "results/bcfs/{sample}.calls.vcf",    
    output:
        callsgz = "results/bcfs/{sample}.calls.vcf.gz",
    log:
        "logs/bgzip/{sample}.log",
    shell:
        """
        bgzip {input.calls} 2> {log}
        """

rule tabix:
    input:
        calls = "results/bcfs/{sample}.calls.vcf.gz",    
    output:
        calls_tbi = "results/bcfs/{sample}.calls.vcf.gz.tbi",
    conda:
        "../envs/AmpSeq.yaml"
    log:
        "logs/tabix/{sample}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """

# wont work with over 1000 files 
rule bcftools_merge:
    input:
        bcfs = expand("results/bcfs/{sample}.calls.vcf.gz", sample=samples),
        idx = expand("results/bcfs/{sample}.calls.vcf.gz.tbi", sample=samples)
    output:
        vcf = "results/vcfs/AgamDaoLSTM_merged.vcf",
    log:
        "logs/bcftools_merge.log",
    conda:
        "../envs/AmpSeq.yaml"
    threads: 12
    shell:
        """
        bcftools merge --threads {threads} -o {output.vcf} -O v {input.bcfs} 2> {log}
        """