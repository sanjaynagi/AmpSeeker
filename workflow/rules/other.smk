# An example collection of Snakemake rules imported in the main Snakefile.

#rule demultiplex:
#    input:
#        read1, read2 = expand("resources/reads/AGAMDAO_{n}.fastq.gz", n=[1,2])
 #   output:
  #      #temp("data/tempalignments/sample.bam")
#    log:
   #     align="logs/align_bwa/{sample}.log",
   ##     sort="logs/sort/{sample}.log",
   #     samblaster="logs/samblaster/{sample}.log.log"
#    resources:bwa=1
#    threads:8
 #   params:
  #      ref="data/reference/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa",
#        tag="'@RG\\tID:{sample}\\tSM:{sample}\\tPU:nye\\tPL:nye\\tLB:{sample}_lb{sample}'"
#    shell:
 #       """
 #       bwa mem -t {threads} {params.ref} {input.reads} -R {params.tag} 2> {log.align} |
 #       samblaster 2> {log.samblaster} | samtools sort -@{threads} -o {output} 2> {log.sort}
 #       """
#
# demuxbyname.sh in=AGAMDAO_LSTM2.R1.fq.gz in2=AGAMDAO_LSTM2.R2.fq.gz out=%_#AGAMDAO.fq.gz outu=AGAMDAO_unmatched.fq.gz barcode=T stats=readStatistics.txt
