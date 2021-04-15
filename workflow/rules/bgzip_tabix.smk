
"""
A snakemake rule file which includes rules for bgzipping, tabix indexing, and merging bcf/vcfs.
Kept in a separate rule file due to maintain visibility of the main, analysis rules in analysis.smk
"""

rule bgzip:
    input:
        calls = "results/bcfs/{sample}.calls.bcf",    
    output:
        callsgz = "results/bcfs/{sample}.calls.bcf.gz",
    log:
        "logs/bgzip/{sample}.log",
    shell:
        """
        bgzip {input.calls} 2> {log}
        """

rule tabix:
    input:
        calls = "results/bcfs/{sample}.calls.bcf.gz",    
    output:
        calls_tbi = "results/bcfs/{sample}.calls.bcf.gz.tbi",
    log:
        "logs/tabix/{sample}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """

# wont work with over 1000 files 
rule bcftools_merge1:
    input:
        bcfs = expand("results/bcfs/{sample}.calls.bcf.gz", sample=samples1),
        idx = expand("results/bcfs/{sample}.calls.bcf.gz.tbi", sample=samples1)
    output:
        vcf = "results/vcfs/AgamDaoLSTM1.vcf",
    log:
        "logs/bcftools_merge.log",
    threads: 12
    shell:
        """
        bcftools merge --threads {threads} -o {output.vcf} -O v {input.bcfs} 2> {log}
        """

rule bcftools_merge2:
    input:
        bcfs = expand("results/bcfs/{sample}.calls.bcf.gz", sample=samples2),
        idx = expand("results/bcfs/{sample}.calls.bcf.gz.tbi", sample=samples2)
    output:
        vcf = "results/vcfs/AgamDaoLSTM2.vcf",
    log:
        "logs/bcftools_merge2.log",
    threads: 12
    shell:
        """
        bcftools merge --threads {threads} -o {output.vcf} -O v {input.bcfs} 2> {log}
        """


rule bgzip2:
    input:
        vcf = "results/vcfs/AgamDaoLSTM{n}.vcf",
    output:
        vcfgz = "results/vcfs/AgamDaoLSTM{n}.vcf.gz",
    log:
        "logs/bgzip/main_{n}.log",
    shell:
        """
        bgzip {input.vcf} 2> {log}
        """

rule tabix2:
    input:
        vcfgz = "results/vcfs/AgamDaoLSTM{n}.vcf.gz",
    output:
        tbi = "results/vcfs/AgamDaoLSTM{n}.vcf.gz.tbi"
    log:
        "logs/tabix/main_{n}.log",
    shell:
        """
        tabix {input.vcfgz} 2> {log}
        """

rule bcftools_merge3:
    input:
        vcf = expand("results/vcfs/AgamDaoLSTM{n}.vcf.gz", n=[1,2]),
        tbi = expand("results/vcfs/AgamDaoLSTM{n}.vcf.gz.tbi", n=[1,2])
    output:
        vcf = "results/vcfs/AgamDaoLSTM_merged.vcf",
    log:
        "logs/bcftools_merge3.log",
    shell:
        """
        bcftools merge -o {output.vcf} -O v {input} 2> {log}
        """