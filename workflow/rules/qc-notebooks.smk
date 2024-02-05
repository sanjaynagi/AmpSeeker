
rule reads_per_well:
    input:
        nb=f"{workflow.basedir}/notebooks/reads-per-well.ipynb",
        kernel="results/.kernel.set",
        bam=expand("results/alignments/{sample}.bam", sample=samples),
        bai=expand("results/alignments/{sample}.bam.bai", sample=samples),
        stats=expand("results/alignments/bamStats/{sample}.flagstat", sample=samples),
        metadata=config["metadata"],
    output:
        nb="results/notebooks/reads-per-well.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/reads-per-well.ipynb",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/reads-per-well.log",
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p metadata_path {input.metadata} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


rule read_qc:
    input:
        nb=f"{workflow.basedir}/notebooks/read-quality.ipynb",
        kernel="results/.kernel.set",
        fastp=expand("results/qc/fastp_reports/{sample}.json", sample=samples),
        metadata=config["metadata"],
    output:
        nb="results/notebooks/read-quality.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/read-quality.ipynb",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/reads-quality.log",
    params:
        index_qc=config["bcl-convert"],
        wd=wkdir,
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p metadata_path {input.metadata} -p index_read_qc {params.index_qc} -p wkdir {params.wd} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


rule coverage:
    input:
        nb=f"{workflow.basedir}/notebooks/coverage.ipynb",
        kernel="results/.kernel.set",
        per_base=expand("results/coverage/{sample}.per-base.bed.gz", sample=samples),
        metadata=config["metadata"],
        targets=config["targets"],
    output:
        nb="results/notebooks/coverage.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/coverage.ipynb",
        excel="results/coverage/amplicon_by_sample_depth.xlsx",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/coverage.log",
    params:
        wkdir=wkdir,
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p metadata_path {input.metadata} -p bed_targets_path {input.targets} -p wkdir {params.wkdir} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """
