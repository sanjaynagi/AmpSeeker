
rule fastp_illumina:
    input:
        sample=lambda w: get_fastqs(wildcards=w, platform=config['platform']),
    output:
        trimmed=[
            "results/trimmed-reads/{sample}_1.fastq.gz",
            "results/trimmed-reads/{sample}_2.fastq.gz",
        ],
        html="results/qc/fastp_reports/{sample}.html",
        json="results/qc/fastp_reports/{sample}.json",
    log:
        "logs/fastp/{sample}.log",
    threads: 4
    wrapper:
        "v1.25.0/bio/fastp"

rule bwa_index:
    input:
        ref=config["reference-fasta"],
    output:
        amb=config["reference-fasta"] + ".amb",
        ann=config["reference-fasta"] + ".ann",
        bwt=config["reference-fasta"] + ".bwt",
        pac=config["reference-fasta"] + ".pac",
        sa=config["reference-fasta"] + ".sa",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    log:
        "logs/bwa_index.log",
    shell:
        """
        bwa index {input.ref} 2> {log}
        """

rule bwa_align:
    """
    Align with bwa mem, and sorting by coordinate with samtools sort. then index with samtools.  
    """
    input:
        reads=expand("results/trimmed-reads/{{sample}}_{n}.fastq.gz", n=[1, 2]),
        ref=config["reference-fasta"],
        amb=config["reference-fasta"] + ".amb",
        ann=config["reference-fasta"] + ".ann",
        bwt=config["reference-fasta"] + ".bwt",
        pac=config["reference-fasta"] + ".pac",
        sa=config["reference-fasta"] + ".sa",
    output:
        bam="results/alignments/{sample}.bam",
    log:
        align="logs/bwa_align/{sample}.log",
        sort="logs/sort/{sample}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    resources:
        bwa=1,
    threads: 1  # each sample tiny so perhaps better to run each on single thread
    params:
        tag="'@RG\\tID:{sample}\\tSM:{sample}\\tPU:nye\\tPL:nye\\tLB:{sample}_lb{sample}'",
    shell:
        """
        bwa mem -t {threads} -R {params.tag} {input.ref} {input.reads} 2> {log.align} |
        samtools sort -@{threads} -o {output} 2> {log.sort}
        """


rule mpileup_call_targets:
    """
    Get pileup of reads at target loci and pipe output to bcftoolsCall
    """
    input:
        bam="results/alignments/{sample}.bam",
        index="results/alignments/{sample}.bam.bai",
        reference=config["reference-fasta"],
    output:
        calls="results/vcfs/targets/{sample}.calls.vcf",
    log:
        mpileup="logs/mpileup/targets/{sample}.log",
        call="logs/bcftools_call/targets/{sample}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    params:
        ref=config["reference-fasta"],
        regions=config["targets"],
        depth=2000,
    shell:
        """
        bcftools mpileup -Ov -f {params.ref} -R {params.regions} -a AD --max-depth {params.depth} {input.bam} 2> {log.mpileup} |
        bcftools call -f GQ,GP -m -Ov 2> {log.call} | bcftools sort -Ov -o {output.calls} 2> {log.call}
        """


rule mpileup_call_amplicons:
    """
    Get pileup of reads at target loci and pipe output to bcftoolsCall
    """
    input:
        bam="results/alignments/{sample}.bam",
        index="results/alignments/{sample}.bam.bai",
        reference=config["reference-fasta"],
    output:
        calls="results/vcfs/amplicons/{sample}.calls.vcf",
    log:
        mpileup="logs/mpileup/amplicons/{sample}.log",
        call="logs/bcftools_call/amplicons/{sample}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    params:
        ref=config["reference-fasta"],
        depth=2000,
    shell:
        """
        bcftools mpileup -Ov -I -f {params.ref} -a AD --max-depth {params.depth} {input.bam} 2> {log.mpileup} |
        bcftools call -m -Ov 2> {log.call} | bcftools sort -Ov -o {output.calls} 2> {log.call}
        """

