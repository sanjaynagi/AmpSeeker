
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




rule create_ploidy_file:
    """
    Create a bcftools ploidy file (see map-call-illumina.smk for rationale).
    Nanopore supports ploidy 1 or 2 only; polyploid is rejected at config
    validation in the Snakefile.
    """
    output:
        ploidy_file="results/config/bcftools_ploidy.txt",
    params:
        ploidy=config["ploidy"],
    log:
        "logs/create_ploidy_file.log",
    run:
        with open(output.ploidy_file, "w") as f:
            f.write(f"* * * M {params.ploidy}\n")
            f.write(f"* * * F {params.ploidy}\n")


rule mpileup_call_targets:
    """
    Get pileup of reads at target loci and pipe output to bcftoolsCall
    """
    input:
        bam="results/alignments/{sample}.bam",
        index="results/alignments/{sample}.bam.bai",
        reference=config["reference-fasta"],
        ploidy_file="results/config/bcftools_ploidy.txt",
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
        bcftools mpileup -X ont-sup -Ov -f {params.ref} -R {params.regions} -a AD --max-depth {params.depth} {input.bam} 2> {log.mpileup} |
        bcftools call -f GQ,GP -m --ploidy-file {input.ploidy_file} -Ov 2> {log.call} | bcftools sort -Ov -o {output.calls} 2> {log.call}
        """


rule clair3_call_amplicons:
    """
    Variant calling with Clair3 across entire amplicon regions (discovery mode).
    Runs in gVCF mode so that reference blocks (with depth) are retained. These are needed
    to distinguish wild-type samples from samples with no data when merging (see
    mask_clair3_gvcf and bcftools_merge).
    For haploid organisms, passes --haploid_sensitive (and disables phasing)
    to clair3 so calls are reported as haploid.
    """
    input:
        bam="results/alignments/{sample}.bam",
        bai="results/alignments/{sample}.bam.bai",
        ref=config["reference-fasta"],
        ref_idx=config["reference-fasta"] + ".fai",
        model_path="resources/models/r1041_e82_400bps_sup_v500",
    output:
        gvcf="results/vcfs/amplicons/{sample}.clair3.g.vcf.gz",
    params:
        outdir="results/clair3_tmp/amplicons/{sample}",
        ploidy_flags="--no_phasing_for_fa --haploid_sensitive" if config["ploidy"] == 1 else "",
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
            --gvcf \
            {params.ploidy_flags} \
            --sample_name={wildcards.sample} &> {log}
    
        # Copy and rename output
        cp {params.outdir}/merge_output.gvcf.gz {output.gvcf}
        
        # Clean up temporary directory
        rm -rf {params.outdir}
        """


rule mask_clair3_gvcf:
    """
    Set genotypes to missing (./.) wherever depth is below `min-genotype-depth`.
    Clair3 labels low-/no-coverage reference blocks as 0/0, so without this those
    positions would be counted as wild type. Applies to reference blocks (FORMAT/MIN_DP)
    and variant records (FORMAT/DP).
    """
    input:
        gvcf="results/vcfs/amplicons/{sample}.clair3.g.vcf.gz",
    output:
        vcf="results/vcfs/amplicons/{sample}.calls.vcf",
    params:
        min_dp=config["min-genotype-depth"],
    conda:
        "../envs/AmpSeeker-cli.yaml"
    log:
        "logs/clair3/amplicons/{sample}.mask.log",
    shell:
        """
        bcftools +setGT {input.gvcf} -Ov -o {output.vcf} -- \
            -t q -n . -i '(FMT/MIN_DP<{params.min_dp}) | (FMT/DP<{params.min_dp})' 2> {log}
        """
