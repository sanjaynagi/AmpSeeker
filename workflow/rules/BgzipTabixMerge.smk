"""
A snakemake rule file which includes rules for bgzipping, tabix indexing, and merging bcf/vcfs.
Kept in a separate rule file due to maintain visibility of the main, analysis rules in analysis.smk
"""

rule bgzip:
    input:
<<<<<<< HEAD
        calls = "results/bcfs/{sample}.calls.vcf",    
    output:
        callsgz = "results/bcfs/{sample}.calls.vcf.gz",
    log:
        "logs/bgzip/{sample}.log",
    shell:
        """
        bgzip {input.calls} 2> {log}
        """
=======
        "results/bcfs/{sample}.calls.vcf",    
    output:
        "results/bcfs/{sample}.calls.vcf.gz",
    log:
        "logs/bgzip/{sample}.log",
    wrapper:
        "v1.17.4-17-g62b55d45/bio/bgzip"
>>>>>>> upstream/main

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

rule bcftools_merge:
    input:
<<<<<<< HEAD
        bcfs = expand("results/bcfs/{sample}.calls.vcf.gz", sample=samples1),
        idx = expand("results/bcfs/{sample}.calls.vcf.gz.tbi", sample=samples1)
    output:
        vcf = "results/vcfs/AgamDaoLSTM1.vcf",
    log:
        "logs/bcftools_merge.log",
=======
        bcfs = expand("results/bcfs/{sample}.calls.vcf.gz", sample=samples),
        idx = expand("results/bcfs/{sample}.calls.vcf.gz.tbi", sample=samples)
    output:
        vcf = "results/vcfs/{dataset}.vcf",
    log:
        "logs/bcftools/merge_{dataset}.log",
>>>>>>> upstream/main
    conda:
        "../envs/AmpSeq.yaml"
    threads: 12
    shell:
        """
        bcftools merge --threads {threads} -o {output.vcf} -O v {input.bcfs} 2> {log}
        """

rule bcftools_merge1:
    input:
<<<<<<< HEAD
        vcfs = expand("results/bcfs/{sample}.calls.vcf.gz", sample=samples2),
        idx = expand("results/bcfs/{sample}.calls.vcf.gz.tbi", sample=samples2)
    output:
        vcf = "results/vcfs/AgamDaoLSTM2.vcf",
    log:
        "logs/bcftools_merge2.log",
=======
        bcfs = expand("results/bcfs/{sample}.calls.vcf.gz", sample=samples1),
        idx = expand("results/bcfs/{sample}.calls.vcf.gz.tbi", sample=samples1)
    output:
        vcf = "results/vcfs/{dataset}.1.vcf",
    log:
        "logs/bcftools/merge1_{dataset}.log",
>>>>>>> upstream/main
    conda:
        "../envs/AmpSeq.yaml"
    threads: 12
    shell:
        """
        bcftools merge --threads {threads} -o {output.vcf} -O v {input.bcfs} 2> {log}
        """

rule bcftools_merge2:
    input:
<<<<<<< HEAD
        vcf = "results/vcfs/AgamDaoLSTM{n}.vcf",
    output:
        vcfgz = "results/vcfs/AgamDaoLSTM{n}.vcf.gz",
    log:
        "logs/bgzip/main_{n}.log",
=======
        vcfs = expand("results/bcfs/{sample}.calls.vcf.gz", sample=samples2),
        idx = expand("results/bcfs/{sample}.calls.vcf.gz.tbi", sample=samples2)
    output:
        vcf = "results/vcfs/{dataset}.2.vcf",
    log:
        "logs/bcftools/merge2_{dataset}.log",
>>>>>>> upstream/main
    conda:
        "../envs/AmpSeq.yaml"
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
<<<<<<< HEAD
        vcfgz = "results/vcfs/AgamDaoLSTM{n}.vcf.gz",
    output:
        tbi = "results/vcfs/AgamDaoLSTM{n}.vcf.gz.tbi"
    log:
        "logs/tabix/main_{n}.log",
=======
        vcfgz = "results/vcfs/{dataset}.{n}.vcf.gz",
    output:
        tbi = "results/vcfs/{dataset}.{n}.vcf.gz.tbi"
    log:
        "logs/tabix/{dataset}.{n}.log",
>>>>>>> upstream/main
    conda:
        "../envs/AmpSeq.yaml"
    shell:
        """
        tabix {input.vcfgz} 2> {log}
        """

rule bcftools_merge3:
    input:
<<<<<<< HEAD
        vcf = expand("results/vcfs/AgamDaoLSTM{n}.vcf.gz", n=[1,2]),
        tbi = expand("results/vcfs/AgamDaoLSTM{n}.vcf.gz.tbi", n=[1,2])
    output:
        vcf = "results/vcfs/AgamDaoLSTM_merged.vcf",
    log:
        "logs/bcftools_merge3.log",
=======
        vcf = expand("results/vcfs/{{dataset}}.{n}.vcf.gz", n=[1,2]),
        tbi = expand("results/vcfs/{{dataset}}.{n}.vcf.gz.tbi", n=[1,2]),
    output:
        vcf = "results/vcfs/{dataset}_merged.vcf",
        holder = touch("results/vcfs/.complete.{dataset}.merge_vcfs")
    log:
        "logs/bcftools/merge3_{dataset}.log",
>>>>>>> upstream/main
    conda:
        "../envs/AmpSeq.yaml"
    shell:
        """
        bcftools merge -o {output.vcf} -Ov {input.vcf} 2> {log}
        """
