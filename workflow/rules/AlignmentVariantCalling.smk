# An example collection of Snakemake rules imported in the main Snakefile.

rule GenomeIndex:
    input:
        ref = lambda wildcards: config['ref'][wildcards.ref]
    output:
        idx = touch("resources/reference/.bwa.index.{ref}")
    shell:
        """
        bwa index {input.ref}
        """

rule alignBWA:
    """
    Align with bwa mem, and sorting by coordinate with samtools sort. then index with samtools.  
    """
    input:
        reads = expand("resources/reads/trimmed/{{sample}}_{n}.fq.gz", n=[1,2]),
        ref = lambda wildcards: config['ref'][wildcards.ref],
        idx = "resources/reference/.bwa.index.{ref}"
    output:
        bam = "resources/{ref}/alignments/{sample}.bam"
    log:
        align="logs/align_bwa/{sample}_{ref}.log",
        sort="logs/sort/{sample}_{ref}.log",
    resources:bwa=1
    threads:1 # each sample tiny so perhaps better to run each on single thread
    params:
        tag="'@RG\\tID:{sample}\\tSM:{sample}\\tPU:nye\\tPL:nye\\tLB:{sample}_lb{sample}'"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.reads} -R {params.tag} 2> {log.align} |
        samtools sort -@{threads} -o {output} 2> {log.sort}
        """

rule indexBams:
     input:
        "resources/{ref}/alignments/{sample}.bam"
     output:
        "resources/{ref}/alignments/{sample}.bam.bai"
     log:
        "logs/index_bams/{sample}_index_{ref}.log"
     shell:
        "samtools index {input} {output} 2> {log}"


rule mpileupAndCall:
    """
    Get pileup of reads at target loci and pipe output to bcftoolsCall
    """
    input:
        bam = "resources/{ref}/alignments/{sample}.bam",
        index = "resources/{ref}/alignments/{sample}.bam.bai"
    output:
        calls = "results/{ref}/bcfs/{sample}.calls.vcf"
    log:
        mpileup = "logs/mpileup/{sample}_{ref}.log",
        call = "logs/bcftools_call/{sample}_{ref}.log"
    params:
        ref = lambda wildcards: config['ref'][wildcards.ref],
        regions = lambda wildcards: config['bed'][wildcards.ref],
        depth = 2000
    shell:
        """
        bcftools mpileup -Ov -f {params.ref} -R {params.regions} --max-depth {params.depth} {input.bam} 2> {log.mpileup} |
        bcftools call -m -Ov 2> {log.call} | bcftools sort -Ov -o {output.calls} 2> {log.call}
        """
