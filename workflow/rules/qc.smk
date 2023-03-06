rule TrimFastqs:
  """
  Trim Fastq files with bbduk - remove adapter sequences and low quality bases
  """
    input:
        read1 = "resources/reads/{sample}_1.fastq.gz",
        read2 = "resources/reads/{sample}_2.fastq.gz"
    output:
        read1 = "resources/reads/trimmed/{sample}_1.fastq.gz",
        read2 = "resources/reads/trimmed/{sample}_2.fastq.gz",
        stats = "resources/reads/trimmed/stats/{sample}.txt"
    log:
        "logs/bbduk/{sample}.log"
    params:
        adaptors = "resources/bbtools_adapters.fa",
        trimq = 10,
        qtrimside="rl",
        ktrimside="l",
    shell:
        """
        bbduk.sh in1={input.read1} in2={input.read2} out1={output.read1} out2={output.read2} \
        qtrim={params.qtrimside} trimq={params.trimq} ref={params.adaptors} ktrim={params.ktrimside} rcomp=t \
        stats={output.stats}
        """


rule targetedCoverage:
  """
  Target per-base coverage with mosdepth
  """
    input:
        bam="results/alignments/{sample}.bam",
        idx="results/alignments/{sample}.bam.bai"
    output:
        "results/coverage/{sample}.per-base.bed.gz"
    log:
        "logs/coverage/{sample}.log"
    threads:4
    conda:
        "../envs/AmpSeq_cli.yaml"
    params:
        prefix="results/coverage/{sample}",
        regions = config['bed']
    shell:
        """
        mosdepth {params.prefix} {input.bam} --by {params.regions} --fast-mode --threads {threads} 2> {log}
        """

rule windowedCoverage:
  """
  300 bp windowed coverage with mosdepth
  """
    input:
        bam="results/alignments/{sample}.bam",
        idx="results/alignments/{sample}.bam.bai"
    output:
        "results/wholegenome/coverage/windowed/{sample}.regions.bed.gz"
    log:
        "logs/coverage/windowed_{sample}_wholegenome.log"
    conda:
        "../envs/AmpSeq_cli.yaml"
    threads:4
    params:
        prefix="results/wholegenome/coverage/windowed/{sample}",
    shell:
        """
        mosdepth {params.prefix} {input.bam} --by 300 -n --fast-mode --threads {threads} 2> {log}
        """

rule BamStats:
  """
  Calculate mapping statistics with samtools flagstat
  """
    input:
        bam = "results/alignments/{sample}.bam",
        idx = "results/alignments/{sample}.bam.bai"
    output:
        stats = "results/alignments/bamStats/{sample}.flagstat"
    conda:
        "../envs/AmpSeq_cli.yaml"
    log:
        "logs/BamStats/{sample}.log"
    shell:
        """
        samtools flagstat {input.bam} > {output} 2> {log}
        """
