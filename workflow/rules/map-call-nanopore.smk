# rule nanoplot_qc:
#     """
#     Quality assessment of nanopore reads with NanoPlot
#     """
#     input:
#         fastq="resources/reads/{sample}.fastq.gz",
#     output:
#         html="results/qc/nanoplot/{sample}/NanoPlot-report.html",
#         stats="results/qc/nanoplot/{sample}/NanoStats.txt",
#     params:
#         outdir="results/qc/nanoplot/{sample}",
#     conda:
#         "../envs/AmpSeeker-nanopore.yaml"
#     log:
#         "logs/nanoplot/{sample}.log",
#     threads: 4
#     shell:
#         """
#         NanoPlot --fastq {input.fastq} --outdir {params.outdir} --threads {threads} 2> {log}
#         """

rule fastp_nanopore:
    input:
        sample=lambda w: get_fastqs(wildcards=w, platform=config['platform']),
    output:
        trimmed="results/trimmed-reads/{sample}.fastq.gz",
        html="results/qc/fastp_reports/{sample}.html",
        json="results/qc/fastp_reports/{sample}.json",
    log:
        "logs/fastplong/{sample}.log",
    conda:
        "../envs/AmpSeeker-nanopore.yaml"
    threads: 4
    shell:
        "fastplong -i {input.sample} -o {output.trimmed} 2> {log}"

rule minimap2_index:
    """
    Index reference genome for minimap2
    """
    input:
        ref=config["reference-fasta"],
    output:
        idx=config["reference-fasta"] + ".mmi",
    conda:
        "../envs/AmpSeeker-nanopore.yaml"
    log:
        "logs/minimap2_index.log",
    threads: 4
    shell:
        """
        minimap2 -t {threads} -d {output.idx} {input.ref} 2> {log}
        """


rule minimap2_align:
    """
    Align nanopore reads with minimap2 and sort with samtools
    """
    input:
        reads="results/trimmed-reads/{sample}.fastq.gz",
        ref=config["reference-fasta"],
        idx=config["reference-fasta"] + ".mmi",
    output:
        bam="results/alignments/{sample}.bam",
    log:
        align="logs/minimap2_align/{sample}.log",
        sort="logs/sort/{sample}.log",
    conda:
        "../envs/AmpSeeker-nanopore.yaml"
    threads: 8
    params:
        rg="'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ONT\\tLB:{sample}_lib'",
    shell:
        """
        minimap2 -t {threads} -ax map-ont -R {params.rg} {input.idx} {input.reads} 2> {log.align} | \
        samtools sort -@ {threads} -o {output.bam} 2> {log.sort}
        """


rule clair3_call_amplicons:
    """
    Variant calling with Clair3 across entire amplicon regions (discovery mode)
    """
    input:
        bam="results/alignments/{sample}.bam",
        bai="results/alignments/{sample}.bam.bai",
        ref=config["reference-fasta"],
        ref_idx=config["reference-fasta"] + ".fai",
    output:
        vcf="results/vcfs/amplicons/{sample}.calls.vcf",
    params:
        model_path=config.get("clair3_model", "/opt/models/ont_guppy5_sup"),
        platform="ont",
        outdir="results/clair3_tmp/amplicons/{sample}",
        min_coverage=config.get("min_coverage", 2),
        amplicon_bed=config.get("amplicon_regions", ""),  # Optional: broader amplicon regions
        sample_name="{sample}",
    conda:
        "../envs/AmpSeeker-nanopore.yaml"
    log:
        "logs/clair3/amplicons/{sample}.log",
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}

        run_clair3.sh \
            --bam_fn={input.bam} \
            --ref_fn={input.ref} \
            --threads={threads} \
            --platform={params.platform} \
            --model_path={params.model_path} \
            --output={params.outdir} \
            --min_coverage={params.min_coverage} \
            --sample_name={params.sample_name} 2>> {log}
    
        # Copy and rename output
        cp {params.outdir}/merge_output.vcf.gz {output.vcf}.tmp.gz
        gunzip {output.vcf}.tmp.gz
        mv {output.vcf}.tmp {output.vcf}
        
        # Clean up temporary directory
        rm -rf {params.outdir}
        """