

rule reference_index:
    input:
        ref = config['reference-fasta']
    output:
        idx = config['reference-fasta'] + ".fai"
    conda:
        "../envs/AmpSeeker-cli.yaml"
    log:
        "logs/reference_index.log"
    shell:
        """
        samtools faidx {input.ref} 2> {log}
        """

rule bwa_index:
    input:
        ref = config['reference-fasta']
    output:
        idx = touch("resources/reference/.bwa.index")
    conda:
        "../envs/AmpSeeker-cli.yaml"
    log:
        "logs/bwa_index.log"
    shell:
        """
        bwa index {input.ref} 2> {log}
        """

rule bwa_align:
    """
    Align with bwa mem, and sorting by coordinate with samtools sort. then index with samtools.  
    """
    input:
        reads = expand("results/trimmed-reads/{{sample}}_{n}.fastq.gz", n=[1,2]),
        ref = config['reference-fasta'],
        idx = "resources/reference/.bwa.index"
    output:
        bam = "results/alignments/{sample}.bam"
    log:
        align="logs/bwa_align/{sample}.log",
        sort="logs/sort/{sample}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    resources:bwa=1
    threads:1 # each sample tiny so perhaps better to run each on single thread
    params:
        tag="'@RG\\tID:{sample}\\tSM:{sample}\\tPU:nye\\tPL:nye\\tLB:{sample}_lb{sample}'"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.reads} -R {params.tag} 2> {log.align} |
        samtools sort -@{threads} -o {output} 2> {log.sort}
        """

rule bam_index:
    input:
        "results/alignments/{sample}.bam"
    output:
        "results/alignments/{sample}.bam.bai"
    conda:
        "../envs/AmpSeeker-cli.yaml"
    log:
        "logs/index_bams/{sample}_index.log"
    shell:
        "samtools index {input} {output} 2> {log}"


rule mpileup_call:
    """
    Get pileup of reads at target loci and pipe output to bcftoolsCall
    """
    input:
        bam = "results/alignments/{sample}.bam",
        index = "results/alignments/{sample}.bam.bai"
    output:
        calls = "results/bcfs/{sample}.calls.vcf"
    log:
        mpileup = "logs/mpileup/{sample}.log",
        call = "logs/bcftools_call/{sample}.log"
    conda:
        "../envs/AmpSeeker-cli.yaml"
    params:
        ref = config['reference-fasta'],
        regions = config['targets'],
        depth = 2000
    shell:
        """
        bcftools mpileup -Ov -f {params.ref} -R {params.regions} --max-depth {params.depth} {input.bam} 2> {log.mpileup} |
        bcftools call -m -Ov 2> {log.call} | bcftools sort -Ov -o {output.calls} 2> {log.call}
        """
