


#rule demultiplex:
#    input:
#        read1, read2 = expand("resources/reads/AGAMDAO_{n}.fastq.gz", n=[1,2])
#    resources:bwa=1
#    threads:8

# demuxbyname.sh in=AGAMDAO_LSTM2.R1.fq.gz in2=AGAMDAO_LSTM2.R2.fq.gz out=%_#AGAMDAO.fq.gz outu=AGAMDAO_unmatched.fq.gz barcode=T stats=readStatistics.txt

# Useful demultiplexing info http://protocols.faircloth-lab.org/en/latest/protocols-computer/sequencing/sequencing-demultiplex-a-run.html

rule TrimFastqs:
  """
  Trim Fastq files with bbduk - remove adapter sequences and low quality bases
  """
    input:
        read1 = "resources/reads/{sample}_1.fq.gz",
        read2 = "resources/reads/{sample}_2.fq.gz"
    output:
        read1 = "resources/reads/trimmed/{sample}_1.fq.gz",
        read2 = "resources/reads/trimmed/{sample}_2.fq.gz",
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
        bam="resources/{ref}/alignments/{sample}.bam",
        idx="resources/{ref}/alignments/{sample}.bam.bai"
    output:
        "results/{ref}/coverage/{sample}.per-base.bed.gz"
    log:
        "logs/coverage/{sample}_{ref}.log"
    threads:4
    params:
        prefix="results/{ref}/coverage/{sample}",
        regions = lambda wildcards: config['bed'][wildcards.ref]
    shell:
        """
        mosdepth {params.prefix} {input.bam} --by {params.regions} --fast-mode --threads {threads} 2> {log}
        """

rule windowedCoverage:
  """
  300 bp windowed coverage with mosdepth
  """
    input:
        bam="resources/wholegenome/alignments/{sample}.bam",
        idx="resources/wholegenome/alignments/{sample}.bam.bai"
    output:
        "results/wholegenome/coverage/windowed/{sample}.regions.bed.gz"
    log:
        "logs/coverage/windowed_{sample}_wholegenome.log"
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
        bam = "resources/{ref}/alignments/{sample}.bam",
        idx = "resources/{ref}/alignments/{sample}.bam.bai"
    output:
        stats = "resources/{ref}/alignments/bamStats/{sample}.flagstat"
    log:
        "logs/BamStats/{sample}_{ref}.log"
    shell:
        """
        samtools flagstat {input.bam} > {output} 2> {log}
        """
