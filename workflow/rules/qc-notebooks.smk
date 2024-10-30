
rule run_information:
    input:
        nb=f"{workflow.basedir}/notebooks/run-information.ipynb",
        kernel="results/.kernel.set",
        metadata=config["metadata"],
        targets=config["targets"],
    output:
        nb="results/notebooks/run-information.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/run-information.ipynb",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/run-information.log",
    params:
        wkdir=wkdir,
        dataset=dataset,
        panel=panel,
        cohort_cols=cohort_cols,
        config_path = workflow.configfiles[0]
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p config_path {params.config_path} -p metadata_path {input.metadata} -p bed_targets_path {input.targets} -p panel {params.panel} -p dataset {params.dataset} -p cohort_cols {params.cohort_cols} -p wkdir {params.wkdir} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """

rule run_statistics:
    input:
        nb=f"{workflow.basedir}/notebooks/run-statistics.ipynb",
        kernel="results/.kernel.set",
        metadata=config["metadata"],
        demultiplex_stats = "results/bcl_output/Reports/Demultiplex_Stats.csv",
    output:
        nb="results/notebooks/run-statistics.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/run-statistics.ipynb",
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/run-statistics.log",
    params:
        wkdir=wkdir,
        cohort_cols=cohort_cols,
        plate_info=plate_info
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python -p plate_info {params.plate_info} -p metadata_path {input.metadata} -p cohort_cols {params.cohort_cols} -p wkdir {params.wkdir} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """

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


rule sample_quality_control:
    input:
        nb=f"{workflow.basedir}/notebooks/sample-quality-control.ipynb",
        kernel="results/.kernel.set",
        region=expand("results/coverage/{sample}.regions.bed.gz", sample=samples),
        vcf=expand("results/vcfs/targets/{dataset}.annot.vcf", dataset=dataset),
        metadata=config["metadata"],
        targets=config["targets"],
    output:
        nb="results/notebooks/sample-quality-control.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/sample-quality-control.ipynb",
        qcpass_metadata = "results/config/metadata.qcpass.tsv"
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/sample-quality-control.log",
    params:
        wkdir=wkdir,
        cohort_cols=cohort_cols,
        sample_threshold = config['quality-control']['sample-total-reads-threshold'],
        panel=panel, 
        dataset=dataset
    shell:
        """
        papermill {input.nb} {output.nb} -k AmpSeq_python - dataset {params.dataset} -p panel {params.panel} -p metadata_path {input.metadata} -p cohort_cols {params.cohort_cols} -p bed_targets_path {input.targets} -p vcf_path {input.vcf} -p wkdir {params.wkdir} -p sample_total_read_threshold {params.sample_threshold} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """
