
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
        "fastplong -i {input.sample} -o {output.trimmed} --html {output.html} --json {output.json} 2> {log}"

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



rule clair3_call_targets:
    """
    Variant calling with Clair3 at specific target sites using BED file
    """
    input:
        bam="results/alignments/{sample}.bam",
        bai="results/alignments/{sample}.bam.bai",
        ref=config["reference-fasta"],
        ref_idx=config["reference-fasta"] + ".fai",
        model_path="resources/models/ont_guppy5",
        bed=config["targets"]
    output:
        vcf="results/vcfs/targets/{sample}.calls.vcf",
    params:
        outdir="results/clair3_tmp/targets/{sample}",
        # ploidy = "--no_phasing_for_fa --haploid_sensitive" if ploidy == 1 else "",
    conda:
        "../envs/AmpSeeker-nanopore.yaml"
    log:
        "logs/clair3/targets/{sample}.log",
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        
        run_clair3.sh \
            --bam_fn={input.bam} \
            --ref_fn={input.ref} \
            --threads={threads} \
            --platform=ont \
            --model_path={input.model_path} \
            --output={params.outdir} \
            --include_all_ctgs \
            --bed_fn={input.bed} \
            --sample_name={wildcards.sample} &> {log}
        
        # Copy and rename output
        cp {params.outdir}/merge_output.vcf.gz {output.vcf}.tmp.gz
        gunzip {output.vcf}.tmp.gz
        mv {output.vcf}.tmp {output.vcf}
        
        # Clean up temporary directory
        rm -rf {params.outdir}
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
        model_path="resources/models/ont_guppy5",
    output:
        vcf="results/vcfs/amplicons/{sample}.calls.vcf",
    params:
        outdir="results/clair3_tmp/amplicons/{sample}",
        # ploidy = "--no_phasing_for_fa --haploid_sensitive" if ploidy == 1 else "",
    conda:
        "../envs/AmpSeeker-nanopore.yaml"
    log:
        "logs/clair3/amplicons/{sample}.log",
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}

        run_clair3.sh \
            --bam_fn={input.bam} \
            --ref_fn={input.ref} \
            --threads={threads} \
            --platform=ont \
            --model_path={input.model_path} \
            --output={params.outdir} \
            --include_all_ctgs \
            --sample_name={wildcards.sample} &> {log}
    
        # Copy and rename output
        cp {params.outdir}/merge_output.vcf.gz {output.vcf}.tmp.gz
        gunzip {output.vcf}.tmp.gz
        mv {output.vcf}.tmp {output.vcf}
        
        # Clean up temporary directory
        rm -rf {params.outdir}
        """